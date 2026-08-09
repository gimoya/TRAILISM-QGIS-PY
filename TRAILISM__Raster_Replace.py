"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterRasterDestination,
    QgsProcessingException,
    QgsPointXY,
    QgsCoordinateTransform,
    QgsProject,
    QgsRasterLayer,
)
from qgis.analysis import (
    QgsRasterCalculator,
    QgsRasterCalculatorEntry,
)
import os
import tempfile
import shutil
import processing


class RasterReplaceAlgorithm(QgsProcessingAlgorithm):
    """
    Replaces DSM raster values with DEM raster values within polygon areas.
    Uses efficient built-in QGIS algorithms (rasterize + calculator).
    DSM (Digital Surface Model) includes vegetation and buildings.
    DEM (Digital Elevation Model) represents the actual terrain surface.
    Optionally levels the ground by using a single DEM value (sampled at polygon centroid)
    for all cells within each polygon.
    """

    DSM_RASTER = 'DSM_RASTER'
    DEM_RASTER = 'DEM_RASTER'
    POLYGONS = 'POLYGONS'
    LEVEL_GROUND = 'LEVEL_GROUND'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return RasterReplaceAlgorithm()

    def name(self):
        return 'raster_replace'

    def displayName(self):
        return 'TRAILISM: Raster Replace DSM with DEM'

    def shortHelpString(self):
        return (
            'TRAILISM ALGORITHM:\n\n'
            'Replaces DSM raster values with DEM raster values within polygon areas.\n'
            'DSM (Digital Surface Model) includes vegetation and buildings.\n'
            'DEM (Digital Elevation Model) represents the actual terrain surface.\n\n'
            'Parameters:\n'
            '- DSM raster: Input DSM raster layer (base surface model with vegetation/buildings)\n'
            '- DEM raster: Input DEM raster layer (terrain values to replace DSM with)\n'
            '- Polygon layer: Defines areas where replacement occurs\n'
            '- Level ground: If checked, uses DEM value at polygon centroid for entire polygon.\n'
            '  If unchecked, uses actual DEM values for each cell within polygons.\n\n'
            'Output: New raster with DSM values replaced by DEM values in polygon areas.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.DSM_RASTER,
                'DSM raster layer (Digital Surface Model)'
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.DEM_RASTER,
                'DEM raster layer (Digital Elevation Model)'
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.POLYGONS,
                'Polygon layer - defines areas for replacement',
                [QgsProcessing.TypeVectorPolygon],
                optional=False
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.LEVEL_GROUND,
                'Level ground - use DEM value at polygon centroid for entire polygon',
                defaultValue=False
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT,
                'Output raster'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        dsm_layer = self.parameterAsRasterLayer(parameters, self.DSM_RASTER, context)
        if dsm_layer is None:
            raise QgsProcessingException("Invalid DSM raster layer")

        dem_layer = self.parameterAsRasterLayer(parameters, self.DEM_RASTER, context)
        if dem_layer is None:
            raise QgsProcessingException("Invalid DEM raster layer")

        polygon_source = self.parameterAsSource(parameters, self.POLYGONS, context)
        if polygon_source is None:
            raise QgsProcessingException("Polygon layer is required")

        level_ground = self.parameterAsBool(parameters, self.LEVEL_GROUND, context)
        output_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)

        # Create temporary directory
        temp_dir = tempfile.mkdtemp()
        try:
            # Step 1: Copy DSM to temp (preserve original)
            dsm_temp = os.path.join(temp_dir, 'dsm_copy.tif')
            feedback.setProgressText("Copying DSM raster...")
            processing.run(
                "gdal:translate",
                {
                    'INPUT': dsm_layer.source(),
                    'TARGET_CRS': None,
                    'NODATA': None,
                    'COPY_SUBDATASETS': False,
                    'OPTIONS': '',
                    'EXTRA': '',
                    'OUTPUT': dsm_temp
                },
                feedback=feedback,
                context=context
            )

            # Step 2: Prepare polygons (reproject if needed)
            dsm_crs = dsm_layer.crs()
            polygon_crs = polygon_source.sourceCrs()
            
            if polygon_crs != dsm_crs:
                polygons_reproj = os.path.join(temp_dir, 'polygons_reproj.gpkg')
                processing.run(
                    "native:reprojectlayer",
                    {
                        'INPUT': parameters[self.POLYGONS],
                        'TARGET_CRS': dsm_crs,
                        'OUTPUT': polygons_reproj
                    },
                    feedback=feedback,
                    context=context
                )
                polygon_layer_path = polygons_reproj
            else:
                # Save polygons to temp file for rasterize
                polygon_layer_path = os.path.join(temp_dir, 'polygons.gpkg')
                processing.run(
                    "native:savefeatures",
                    {
                        'INPUT': parameters[self.POLYGONS],
                        'OUTPUT': polygon_layer_path
                    },
                    feedback=feedback,
                    context=context
                )

            # Step 3: Handle level_ground option
            if level_ground:
                feedback.pushWarning("Level ground mode: using first polygon's centroid value for all polygons")
                dem_dp = dem_layer.dataProvider()
                dem_to_use = dem_layer.source()
                
                # Get first polygon's centroid and create modified DEM
                for feat in polygon_source.getFeatures():
                    geom = feat.geometry()
                    if geom and not geom.isEmpty():
                        if polygon_crs != dsm_crs:
                            xform = QgsCoordinateTransform(polygon_crs, dsm_crs, QgsProject.instance())
                            geom.transform(xform)
                        centroid = geom.centroid()
                        point = centroid.asPoint()
                        dem_val, ok = dem_dp.sample(QgsPointXY(point.x(), point.y()), 1)
                        if ok and dem_val is not None:
                            centroid_value = float(dem_val)
                            
                            # Create constant raster with centroid value
                            constant_raster = os.path.join(temp_dir, 'constant_dem.tif')
                            processing.run(
                                "gdal:rasterize",
                                {
                                    'INPUT': polygon_layer_path,
                                    'BURN': centroid_value,
                                    'UNITS': 0,
                                    'WIDTH': dsm_layer.width(),
                                    'HEIGHT': dsm_layer.height(),
                                    'EXTENT': dsm_layer.extent(),
                                    'NODATA': None,
                                    'OPTIONS': '',
                                    'DATA_TYPE': 5,
                                    'INIT': None,
                                    'INVERT': False,
                                    'EXTRA': '',
                                    'OUTPUT': constant_raster
                                },
                                feedback=feedback,
                                context=context
                            )
                            
                            # Create mask
                            mask_raster = os.path.join(temp_dir, 'mask.tif')
                            processing.run(
                                "gdal:rasterize",
                                {
                                    'INPUT': polygon_layer_path,
                                    'BURN': 1,
                                    'UNITS': 0,
                                    'WIDTH': dsm_layer.width(),
                                    'HEIGHT': dsm_layer.height(),
                                    'EXTENT': dsm_layer.extent(),
                                    'NODATA': 0,
                                    'OPTIONS': '',
                                    'DATA_TYPE': 1,
                                    'INIT': 0,
                                    'INVERT': False,
                                    'EXTRA': '',
                                    'OUTPUT': mask_raster
                                },
                                feedback=feedback,
                                context=context
                            )
                            
                            # Merge: (mask == 1) * constant + (mask == 0) * DEM
                            dem_to_use = os.path.join(temp_dir, 'dem_modified.tif')
                            processing.run(
                                "qgis:rastercalculator",
                                {
                                    'EXPRESSION': '((A@1 == 1) * B@1) + ((A@1 == 0) * C@1)',
                                    'LAYERS': [mask_raster, constant_raster, dem_layer.source()],
                                    'CELLSIZE': 0,
                                    'EXTENT': None,
                                    'CRS': None,
                                    'OUTPUT': dem_to_use
                                },
                                feedback=feedback,
                                context=context
                            )
                            break
            else:
                dem_to_use = dem_layer.source()

            # Step 4: Rasterize polygons with -9999 into DSM copy
            feedback.setProgressText("Rasterizing polygon mask into DSM...")
            processing.run(
                "gdal:rasterize_over_fixed_value",
                {
                    'INPUT': polygon_layer_path,
                    'INPUT_RASTER': dsm_temp,
                    'BURN': -9999,
                    'ADD': False,
                    'EXTRA': ''
                },
                feedback=feedback,
                context=context
            )

            # Step 5: Merge with raster calculator
            # Formula: (DSM = -9999) ? DEM : DSM
            feedback.setProgressText("Merging rasters with calculator...")
            
            # Load rasters as layers
            dsm_calc_layer = QgsRasterLayer(dsm_temp, 'dsm_temp')
            dem_calc_layer = QgsRasterLayer(dem_to_use, 'dem_temp')
            
            if not dsm_calc_layer.isValid():
                raise QgsProcessingException(f"Cannot load DSM raster for calculator: {dsm_temp}")
            if not dem_calc_layer.isValid():
                raise QgsProcessingException(f"Cannot load DEM raster for calculator: {dem_to_use}")

            # Create calculator entries
            entries = []
            
            dsm_entry = QgsRasterCalculatorEntry()
            dsm_entry.ref = 'A@1'
            dsm_entry.raster = dsm_calc_layer
            dsm_entry.bandNumber = 1
            entries.append(dsm_entry)
            
            dem_entry = QgsRasterCalculatorEntry()
            dem_entry.ref = 'B@1'
            dem_entry.raster = dem_calc_layer
            dem_entry.bandNumber = 1
            entries.append(dem_entry)
            
            # Create calculator with transformContext (required since QGIS 3.8)
            formula = '(A@1 = -9999) * B@1 + (A@1 != -9999) * A@1'
            calc = QgsRasterCalculator(
                formula,
                output_path,
                'GTiff',
                dsm_calc_layer.extent(),
                dsm_calc_layer.width(),
                dsm_calc_layer.height(),
                entries,
                context.transformContext()
            )
            
            # Process calculation with feedback support
            result = calc.processCalculation(feedback)
            if result != QgsRasterCalculator.Success:
                error_msg = calc.lastError() if hasattr(calc, 'lastError') else f"Calculator error code: {result}"
                raise QgsProcessingException(f"Raster calculator failed: {error_msg}")

            # Copy style from original DSM to output
            feedback.setProgressText("Copying style from DSM raster...")
            output_layer = QgsRasterLayer(output_path, 'output')
            if output_layer.isValid() and dsm_layer.renderer():
                try:
                    # Copy renderer from DSM to output
                    output_layer.setRenderer(dsm_layer.renderer().clone())
                    output_layer.triggerRepaint()
                    
                    # Also copy full style using QgsMapLayerStyle
                    style = QgsMapLayerStyle()
                    style.readFromLayer(dsm_layer)
                    style.writeToLayer(output_layer)
                except Exception:
                    pass

            feedback.setProgress(100)
            return {self.OUTPUT: output_path}

        finally:
            try:
                shutil.rmtree(temp_dir)
            except Exception:
                pass

