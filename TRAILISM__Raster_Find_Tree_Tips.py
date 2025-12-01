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
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterExtent,
    QgsProcessingException,
    QgsPointXY,
    QgsRasterFileWriter,
    QgsRasterBlock,
    QgsRasterPipe,
    QgsRasterLayer,
    QgsSingleBandPseudoColorRenderer,
    QgsColorRampShader,
    QgsRasterShader,
    QgsRasterBandStats,
    Qgis,
)
from qgis.PyQt.QtGui import QColor
import math
import os


class RasterFindTreeTipsAlgorithm(QgsProcessingAlgorithm):
    """
    Creates a minimum slope raster map: for each cell, calculates average slope
    in 8 directions and outputs the minimum slope value. Positive values indicate
    peaks (downhill in all directions), negative values indicate depressions or slopes.
    """

    INPUT_RASTER = 'INPUT_RASTER'
    EXTENT = 'EXTENT'
    LOOKING_DISTANCE = 'LOOKING_DISTANCE'
    MIN_THRESHOLD = 'MIN_THRESHOLD'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return RasterFindTreeTipsAlgorithm()

    def name(self):
        return 'raster_find_tree_tips'

    def displayName(self):
        return 'TRAILISM: Raster Find Tree Tips'

    def shortHelpString(self):
        return (
            'TRAILISM ALGORITHM:\n\n'
            'Creates an average slope angle raster map for tree canopy detection.\n'
            'For each cell, calculates average slope angle in 8 directions by stepping\n'
            'through all cells up to the looking distance. Outputs the average\n'
            'slope angle in degrees if above minimum threshold. Positive values = peaks.\n\n'
            'Parameters:\n'
            '- Extent: select on canvas or use current canvas view (optional)\n'
            '- Looking distance: distance to check in each direction (meters).\n'
            '  Should be 1.5-2x typical crown diameter to ensure comparison to ground.\n'
            '  For coniferous trees (5-12m crowns): use 10-20m (default: 12m)\n'
            '- Minimum threshold: minimum average slope angle in degrees.\n'
            '  Typical values: 10-15° (small trees), 15-25° (medium), 25-35° (large)\n'
            '  Default: 15° (filters small peaks, detects medium+ trees)\n\n'
            'Output: raster layer with average slope angles in degrees above threshold.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_RASTER,
                'Input raster layer'
            )
        )

        self.addParameter(
            QgsProcessingParameterExtent(
                self.EXTENT,
                'Extent - select on canvas or use current canvas view',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.LOOKING_DISTANCE,
                'Looking distance (meters) - should be 1.5-2x typical crown diameter',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=12.0,
                minValue=0.001
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MIN_THRESHOLD,
                'Minimum threshold (degrees) - typical: 10-25° for coniferous trees',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=15.0,
                minValue=0.0,
                maxValue=90.0
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT,
                'Output slope raster'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Main processing algorithm: calculates average slope angle raster map.
        
        For each cell, calculates average slope angle in 8 directions (N, NE, E, SE, S, SW, W, NW)
        over a specified looking distance. Outputs the average slope angle in degrees
        if above minimum threshold, otherwise sets to nodata.
        Positive values indicate peaks, negative values indicate depressions or slopes.
        """
        # Get input raster layer
        raster_layer = self.parameterAsRasterLayer(parameters, self.INPUT_RASTER, context)
        if raster_layer is None:
            raise QgsProcessingException("Invalid raster layer")

        # Get looking distance parameter (distance to check in each direction, in meters)
        looking_distance = self.parameterAsDouble(parameters, self.LOOKING_DISTANCE, context)
        if looking_distance <= 0:
            raise QgsProcessingException("Looking distance must be positive")

        # Get minimum threshold parameter (degrees)
        min_threshold = self.parameterAsDouble(parameters, self.MIN_THRESHOLD, context)
        if min_threshold < 0 or min_threshold > 90:
            raise QgsProcessingException("Minimum threshold must be between 0 and 90 degrees")

        # Get raster properties
        dp = raster_layer.dataProvider()
        raster_extent = raster_layer.extent()
        crs = raster_layer.crs()
        
        # Try multiple methods to get nodata value
        nodata = None
        if hasattr(dp, "sourceNoDataValue"):
            try:
                nodata = dp.sourceNoDataValue(1)
            except:
                pass
        if nodata is None and hasattr(dp, "userNoDataValue"):
            try:
                nodata = dp.userNoDataValue(1)
            except:
                pass
        if nodata is None and hasattr(raster_layer, "dataProvider"):
            try:
                nodata = raster_layer.dataProvider().sourceNoDataValue(1)
            except:
                pass
        
        feedback.pushInfo(f"Raster nodata value: {nodata if nodata is not None else 'not set'}")

        # Get extent (optional) - automatically transformed to raster CRS
        # If provided, only cells within this extent will be processed
        extent_rect = self.parameterAsExtent(parameters, self.EXTENT, context, crs)
        
        if extent_rect and not extent_rect.isEmpty():
            feedback.pushInfo(f"Extent filter: {extent_rect.xMinimum():.2f}, {extent_rect.yMinimum():.2f} to {extent_rect.xMaximum():.2f}, {extent_rect.yMaximum():.2f}")
            feedback.pushInfo(f"Raster extent: {raster_extent.xMinimum():.2f}, {raster_extent.yMinimum():.2f} to {raster_extent.xMaximum():.2f}, {raster_extent.yMaximum():.2f}")
            # Check if extent overlaps with raster
            if not raster_extent.intersects(extent_rect):
                raise QgsProcessingException(f"Provided extent does not overlap with raster extent. Raster: {raster_extent.xMinimum():.2f}, {raster_extent.yMinimum():.2f} to {raster_extent.xMaximum():.2f}, {raster_extent.yMaximum():.2f}")

        # Get raster dimensions
        width = raster_layer.width()
        height = raster_layer.height()
        
        if width <= 0 or height <= 0:
            raise QgsProcessingException("Invalid raster dimensions")

        # Get cell size in map units (meters)
        cell_size_x_map = abs(raster_layer.rasterUnitsPerPixelX())
        cell_size_y_map = abs(raster_layer.rasterUnitsPerPixelY())
        avg_cell_size = (cell_size_x_map + cell_size_y_map) / 2.0

        # Convert looking distance from meters to number of cells
        # This determines how many cells to step through in each direction
        distance_cells = int(round(looking_distance / avg_cell_size))
        if distance_cells < 1:
            distance_cells = 1

        # Pre-calculate border limits: cells within distance_cells of edges cannot be processed
        # because we can't check all 8 directions for them
        min_col = distance_cells
        max_col = width - distance_cells
        min_row = distance_cells
        max_row = height - distance_cells
        
        # Pre-calculate coordinate conversion constants (for extent filtering if needed)
        xmin = raster_extent.xMinimum()
        ymax = raster_extent.yMaximum()

        # Read entire raster
        feedback.setProgressText(f"Reading raster data ({width} × {height} cells)...")
        block = dp.block(1, raster_extent, width, height)
        
        if block is None or block.isEmpty():
            raise QgsProcessingException("Failed to read raster data")

        # Extract raster values into dictionary for efficient access
        # Include ALL cells (including borders) - border cells are needed for neighbor lookups
        # when calculating slopes for non-border cells, even though border cells themselves
        # won't be processed (they'll be set to nodata in output)
        raster_data = {}  # (col, row) -> value
        cells_checked = 0
        cells_valid = 0
        cells_in_extent = 0
        cells_nodata = 0
        cells_nan = 0
        sample_values = []  # Collect first few values for debugging
        
        for row in range(height):
            if feedback.isCanceled():
                return {self.OUTPUT: None}
            for col in range(width):
                cells_checked += 1
                try:
                    val = block.value(row, col)
                    if val is not None:
                        try:
                            zv = float(val)
                            
                            # Check for NaN
                            if math.isnan(zv):
                                cells_nan += 1
                                continue
                            
                            # Check for infinity
                            if math.isinf(zv):
                                continue
                            
                            # Check if value is valid (not nodata)
                            is_nodata = False
                            if nodata is not None:
                                # Use relative tolerance for nodata comparison
                                if abs(zv - nodata) < max(1e-6, abs(nodata) * 1e-6):
                                    is_nodata = True
                            
                            if not is_nodata:
                                cells_valid += 1
                                
                                # Collect sample values for debugging (first 10)
                                if len(sample_values) < 10:
                                    sample_values.append(zv)
                                
                                # Only convert coordinates if extent filtering is needed
                                if extent_rect and not extent_rect.isEmpty():
                                    x = xmin + (col + 0.5) * cell_size_x_map
                                    y = ymax - (row + 0.5) * cell_size_y_map
                                    if not extent_rect.contains(QgsPointXY(x, y)):
                                        continue
                                    cells_in_extent += 1
                                
                                raster_data[(col, row)] = zv
                            else:
                                cells_nodata += 1
                        except (ValueError, TypeError) as e:
                            pass
                except Exception:
                    pass

        feedback.pushInfo(f"Raster reading: {cells_checked} cells checked")
        feedback.pushInfo(f"  - Valid cells: {cells_valid}")
        feedback.pushInfo(f"  - Nodata cells: {cells_nodata}")
        feedback.pushInfo(f"  - NaN cells: {cells_nan}")
        if sample_values:
            feedback.pushInfo(f"  - Sample values: {sample_values[:5]}")
        if extent_rect and not extent_rect.isEmpty():
            feedback.pushInfo(f"Extent filtering: {cells_in_extent} cells within extent")

        if not raster_data:
            error_msg = f"No valid raster cells found. "
            if cells_valid == 0:
                if cells_nodata > 0:
                    error_msg += f"All {cells_nodata} cells are marked as nodata. "
                if cells_nan > 0:
                    error_msg += f"Found {cells_nan} NaN values. "
                error_msg += "Check raster nodata settings or raster data validity."
            elif extent_rect and not extent_rect.isEmpty():
                error_msg += f"Found {cells_valid} valid cells in raster, but none within extent bounds."
            raise QgsProcessingException(error_msg)

        feedback.pushInfo(f"Calculating slope angles: {looking_distance:.2f}m distance ({distance_cells} cells), min threshold: {min_threshold:.2f}°")

        # 8 direction vectors (unit vectors for stepping through cells)
        # These define the 8 cardinal and diagonal directions
        direction_vectors = [
            (0, -1),    # N  (north)
            (1, -1),    # NE (northeast)
            (1, 0),     # E  (east)
            (1, 1),     # SE (southeast)
            (0, 1),     # S  (south)
            (-1, 1),    # SW (southwest)
            (-1, 0),    # W  (west)
            (-1, -1)    # NW (northwest)
        ]

        # Pre-calculate Euclidean distances for each direction (constant for all cells)
        # For diagonal directions, the actual distance is longer than the cell count
        # This ensures accurate slope calculations regardless of direction
        euclidean_distances = []
        for dx, dy in direction_vectors:
            dist_x = dx * distance_cells * cell_size_x_map
            dist_y = dy * distance_cells * cell_size_y_map
            euclidean_dist = math.sqrt(dist_x**2 + dist_y**2)
            euclidean_distances.append(euclidean_dist)

        # Create output raster block and initialize with nodata
        output_block = QgsRasterBlock(Qgis.DataType.Float32, width, height)
        output_nodata = -9999.0
        output_block.setNoDataValue(output_nodata)
        
        # Initialize all cells to nodata
        for row in range(height):
            for col in range(width):
                output_block.setValue(row, col, output_nodata)

        total_cells = len(raster_data)
        processed = 0
        progress_interval = max(1000, total_cells // 100)

        feedback.setProgressText("Calculating average slope angles for each cell...")

        # Calculate minimum slope for each cell
        for (col, row), cell_elev in raster_data.items():
            if feedback.isCanceled():
                break

            processed += 1
            if processed % progress_interval == 0:
                progress = int(100 * processed / total_cells)
                feedback.setProgress(progress)

            # Skip border cells for processing (set to nodata in output)
            # But border cells are still in raster_data for neighbor lookups
            if col < min_col or col >= max_col or row < min_row or row >= max_row:
                output_block.setValue(row, col, output_nodata)
                continue

            slopes = []
            all_directions_valid = True

            # Calculate average slope for each of 8 directions
            for dir_idx, (dx, dy) in enumerate(direction_vectors):
                path_elevations = []

                # Step through cells in this direction up to distance_cells
                # Collect all elevation values along the path
                for step in range(1, distance_cells + 1):
                    check_col = col + dx * step
                    check_row = row + dy * step

                    # Check bounds - if we hit border, path is incomplete
                    # This should not happen for non-border cells, but check for safety
                    if check_col < 0 or check_col >= width or check_row < 0 or check_row >= height:
                        all_directions_valid = False
                        break

                    # Get elevation (exclude nodata from calculation but continue path)
                    # Nodata cells are skipped in the average but don't invalidate the path
                    check_elev = raster_data.get((check_col, check_row))
                    if check_elev is not None:
                        path_elevations.append(check_elev)

                # If path is incomplete (hit border), mark as invalid
                if not all_directions_valid:
                    break

                # If no valid elevations found in path, invalid
                if len(path_elevations) == 0:
                    all_directions_valid = False
                    break

                # Calculate slope: (start_elev - avg_path_elev) / distance
                # Slope is positive if cell is higher than path average (downhill)
                # Slope is negative if cell is lower than path average (uphill)
                # Use pre-calculated Euclidean distance for accurate diagonal calculations
                euclidean_distance = euclidean_distances[dir_idx]
                avg_path_elev = sum(path_elevations) / len(path_elevations)
                slope = (cell_elev - avg_path_elev) / euclidean_distance
                # Convert slope to angle in degrees
                angle_degrees = math.degrees(math.atan(slope))
                slopes.append(angle_degrees)

            # Output average slope angle in degrees (or nodata if not all directions valid or below/equal threshold)
            if all_directions_valid and len(slopes) == 8:
                avg_angle = sum(slopes) / len(slopes)
                # Only output if strictly above minimum threshold, otherwise set to nodata
                if avg_angle > min_threshold:
                    output_block.setValue(row, col, avg_angle)
                else:
                    output_block.setValue(row, col, output_nodata)
            else:
                # Invalid directions - set to nodata
                output_block.setValue(row, col, output_nodata)

        # Write output raster
        feedback.setProgressText("Writing output raster...")
        output_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)
        
        # Write raster using GDAL directly to ensure proper data type and statistics
        # This approach creates a new raster with Float32 type and computes statistics
        # so QGIS displays the correct data range (slope percentages) instead of original elevation range
        try:
            from osgeo import gdal, osr
            
            # Create new GeoTIFF with Float32 data type
            driver = gdal.GetDriverByName('GTiff')
            if driver is None:
                raise QgsProcessingException("Cannot get GTiff driver")
            
            # Create dataset with Float32 type (single band)
            ds = driver.Create(output_path, width, height, 1, gdal.GDT_Float32)
            if ds is None:
                raise QgsProcessingException("Cannot create output raster file")
            
            # Set geotransform (defines pixel size and location)
            geotransform = [
                raster_extent.xMinimum(),  # top-left x coordinate
                cell_size_x_map,           # pixel width
                0,                         # rotation (0 = north up)
                raster_extent.yMaximum(),  # top-left y coordinate
                0,                         # rotation (0 = north up)
                -cell_size_y_map           # pixel height (negative because y decreases downward)
            ]
            ds.SetGeoTransform(geotransform)
            
            # Set coordinate reference system
            srs = osr.SpatialReference()
            srs.ImportFromWkt(crs.toWkt())
            ds.SetProjection(srs.ExportToWkt())
            
            # Get band and set no data value
            band = ds.GetRasterBand(1)
            band.SetNoDataValue(output_nodata)
            
            # Prepare data array from output_block (convert to flat array)
            import array
            data_array = array.array('f')
            for row in range(height):
                for col in range(width):
                    val = output_block.value(row, col)
                    data_array.append(val if val != output_nodata else output_nodata)
            
            # Write entire block at once (more efficient than row-by-row)
            band.WriteRaster(0, 0, width, height, data_array.tobytes(), width, height, gdal.GDT_Float32)
            band.FlushCache()
            
            # Calculate statistics for proper display in QGIS
            # This ensures QGIS shows the correct min/max range for slope percentages
            band.ComputeStatistics(False)
            
            # Close dataset (ensures file is written and closed properly)
            ds = None
            
        except ImportError:
            # Fallback to QGIS method if GDAL not available
            feedback.pushInfo("GDAL not available, using QGIS raster writer (may show incorrect data range)")
            
            writer = QgsRasterFileWriter(output_path)
            writer.setOutputProviderKey('gdal')
            writer.setOutputFormat('GTiff')
            
            pipe = QgsRasterPipe()
            provider = dp.clone()
            if not pipe.set(provider):
                raise QgsProcessingException("Cannot set pipe provider")
            
            error = writer.writeRaster(pipe, width, height, raster_extent, crs, context.transformContext())
            if error != QgsRasterFileWriter.NoError:
                raise QgsProcessingException(f"Error writing raster: {error}")
            
            temp_layer = QgsRasterLayer(output_path, "temp")
            if not temp_layer.isValid():
                raise QgsProcessingException("Cannot open output raster file")
            
            output_provider = temp_layer.dataProvider()
            if output_provider is None:
                raise QgsProcessingException("Cannot get output raster provider")
            
            output_provider.setEditable(True)
            output_provider.writeBlock(output_block, 1, 0, 0)
            output_provider.setEditable(False)
            output_provider.setNoDataValue(1, output_nodata)

        feedback.pushInfo(f"Created average slope angle raster with {total_cells} cells processed (threshold: {min_threshold:.2f}°)")

        # Style output raster: positive values = red, zero and negative = transparent
        feedback.setProgressText("Applying styling to output raster...")
        output_layer = QgsRasterLayer(output_path, "Slope Prominence")
        if output_layer.isValid():
            # Get actual statistics from the raster to set proper value range
            stats = output_layer.dataProvider().bandStatistics(1, QgsRasterBandStats.All, output_layer.extent(), 0)
            min_val = stats.minimumValue if stats.minimumValue > output_nodata else min_threshold
            max_val = stats.maximumValue if stats.maximumValue > output_nodata else max(min_threshold + 10.0, 45.0)
            
            # Ensure we have a reasonable range
            if max_val <= min_val:
                max_val = min_val + 30.0  # Default range of 30 degrees
            
            feedback.pushInfo(f"Styling raster with value range: {min_val:.2f} to {max_val:.2f} degrees")
            
            # Create color ramp shader
            shader = QgsRasterShader()
            color_ramp = QgsColorRampShader()
            color_ramp.setColorRampType(QgsColorRampShader.Type.Interpolated)
            
            # Create color ramp entries based on actual value range
            # Entry 1: Below threshold -> transparent
            entry1 = QgsColorRampShader.ColorRampItem(min_val - 1.0, QColor(255, 0, 0, 0), "Below threshold")
            
            # Entry 2: At threshold -> transparent
            entry2 = QgsColorRampShader.ColorRampItem(min_threshold, QColor(255, 0, 0, 0), "Threshold")
            
            # Entry 3: Just above threshold -> light red
            entry3 = QgsColorRampShader.ColorRampItem(min_threshold + 0.01, QColor(255, 100, 100, 200), "Low prominence")
            
            # Entry 4: Mid range -> medium red
            mid_val = (min_threshold + max_val) / 2.0
            entry4 = QgsColorRampShader.ColorRampItem(mid_val, QColor(255, 50, 50, 255), "Medium prominence")
            
            # Entry 5: Max value -> solid red
            entry5 = QgsColorRampShader.ColorRampItem(max_val, QColor(255, 0, 0, 255), "High prominence")
            
            color_ramp.setColorRampItemList([entry1, entry2, entry3, entry4, entry5])
            shader.setRasterShaderFunction(color_ramp)
            
            # Create renderer
            renderer = QgsSingleBandPseudoColorRenderer(
                output_layer.dataProvider(),
                1,  # band number
                shader
            )
            output_layer.setRenderer(renderer)
            output_layer.triggerRepaint()
            
            # Add layer to project if we have access to it
            try:
                from qgis.core import QgsProject
                QgsProject.instance().addMapLayer(output_layer)
            except:
                pass  # If we can't add to project, that's okay

        return {self.OUTPUT: output_path}

