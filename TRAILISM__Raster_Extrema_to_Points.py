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

from qgis.PyQt.QtCore import QVariant, QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingException,
    QgsFields,
    QgsField,
    QgsFeature,
    QgsFeatureSink,
    QgsGeometry,
    QgsPointXY,
    QgsWkbTypes,
    QgsMarkerSymbol,
    QgsSingleSymbolRenderer,
    QgsRuleBasedRenderer,
    QgsRectangle,
    QgsCoordinateTransform,
    QgsProject,
)


class RasterExtremaAlgorithm(QgsProcessingAlgorithm):
    """
    Finds minimum and maximum raster values within each polygon feature.
    Outputs points at the locations of these extrema.
    """

    INPUT_RASTER = 'INPUT_RASTER'
    POLYGONS = 'POLYGONS'
    EXTREMA_TYPE = 'EXTREMA_TYPE'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return RasterExtremaAlgorithm()

    def name(self):
        return 'raster_extrema'

    def displayName(self):
        return 'TRAILISM: Raster Extrema to Points'

    def shortHelpString(self):
        return (
            'TRAILISM ALGORITHM:\n\n'
            'Finds minimum and maximum raster values within each polygon feature.\n'
            'Outputs points at the locations of these extrema.\n\n'
            'Parameters:\n'
            '- Extrema type: find maxima (peaks), minima (depressions), or both\n\n'
            'Polygon layer is required - finds extrema within each polygon feature.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_RASTER,
                'Input raster layer'
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.POLYGONS,
                'Polygon layer - restrict sampling to areas within polygons',
                [QgsProcessing.TypeVectorPolygon],
                optional=False
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.EXTREMA_TYPE,
                'Extrema type',
                options=['Maxima only (peaks)', 'Minima only (depressions)', 'Both'],
                allowMultiple=False,
                defaultValue=0
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                'Output points'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        raster_layer = self.parameterAsRasterLayer(parameters, self.INPUT_RASTER, context)
        if raster_layer is None:
            raise QgsProcessingException("Invalid raster layer")

        # Get polygon layer - required
        polygon_source = self.parameterAsSource(parameters, self.POLYGONS, context)
        if polygon_source is None:
            raise QgsProcessingException("Polygon layer is required")
        
        polygon_features = []
        for feat in polygon_source.getFeatures():
            geom = feat.geometry()
            if geom and not geom.isEmpty():
                polygon_features.append((geom, feat))
        
        if not polygon_features:
            raise QgsProcessingException("No valid polygon features found")
        
        extrema_type = self.parameterAsEnum(parameters, self.EXTREMA_TYPE, context)

        find_maxima = extrema_type in [0, 2]  # Maxima only or Both
        find_minima = extrema_type in [1, 2]  # Minima only or Both

        # Get raster properties
        dp = raster_layer.dataProvider()
        raster_extent = raster_layer.extent()
        crs = raster_layer.crs()
        
        # Initialize output sink early for cancellation handling
        fields = QgsFields()
        fields.append(QgsField("extrema_type", QVariant.String))
        fields.append(QgsField("value", QVariant.Double))
        fields.append(QgsField("polygon_id", QVariant.Int))

        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT, context, fields, QgsWkbTypes.Point, crs
        )

        if sink is None:
            raise QgsProcessingException("Invalid output sink")
        nodata = dp.sourceNoDataValue(1) if hasattr(dp, "sourceNoDataValue") else None

        # Get raster dimensions
        width = raster_layer.width()
        height = raster_layer.height()
        
        if width <= 0 or height <= 0:
            raise QgsProcessingException("Invalid raster dimensions")

        cell_size_x_map = abs(raster_layer.rasterUnitsPerPixelX())
        cell_size_y_map = abs(raster_layer.rasterUnitsPerPixelY())

        feedback.pushInfo(f"Raster cell size: {cell_size_x_map:.3f} m × {cell_size_y_map:.3f} m")

        # Get raster origin
        xmin = raster_extent.xMinimum()
        ymax = raster_extent.yMaximum()

        # Transform polygons to raster CRS if needed
        polygon_crs = polygon_source.sourceCrs()
        if polygon_crs != crs:
            xform = QgsCoordinateTransform(polygon_crs, crs, QgsProject.instance())
            transformed_polygons = []
            for poly_geom, feat in polygon_features:
                try:
                    poly_geom.transform(xform)
                    transformed_polygons.append((poly_geom, feat))
                except Exception:
                    pass
            polygon_features = transformed_polygons

        # Convert map coordinates to raster cell indices
        def map_to_raster_cell(x, y):
            """Convert map coordinates to raster cell (col, row)."""
            col = int((x - xmin) / cell_size_x_map)
            row = int((ymax - y) / cell_size_y_map)
            col = max(0, min(width - 1, col))
            row = max(0, min(height - 1, row))
            return col, row

        # Convert raster row/col to map coordinates (cell centroid)
        def cell_centroid_to_map(col, row):
            """Convert raster cell (col, row) to map coordinates of cell centroid."""
            x = xmin + (col + 0.5) * cell_size_x_map
            y = ymax - (row + 0.5) * cell_size_y_map
            return x, y

        # Calculate bounding box of all polygons
        all_poly_bounds = [poly_geom.boundingBox() for poly_geom, _ in polygon_features]
        min_x = min(b.xMinimum() for b in all_poly_bounds)
        max_x = max(b.xMaximum() for b in all_poly_bounds)
        min_y = min(b.yMinimum() for b in all_poly_bounds)
        max_y = max(b.yMaximum() for b in all_poly_bounds)
        
        # Convert to raster cell indices
        col_start, row_top = map_to_raster_cell(min_x, max_y)
        col_end, row_bottom = map_to_raster_cell(max_x, min_y)
        
        # Clamp to raster bounds with buffer
        col_start = max(0, col_start - 1)
        col_end = min(width, col_end + 2)
        row_top = max(0, row_top - 1)
        row_bottom = min(height, row_bottom + 2)
        
        read_col_start = col_start
        read_row_start = row_top
        read_width = col_end - col_start
        read_height = row_bottom - row_top

        # Read raster block covering polygon bounds
        read_extent_rect = QgsRectangle(
            xmin + col_start * cell_size_x_map,
            ymax - row_bottom * cell_size_y_map,
            xmin + col_end * cell_size_x_map,
            ymax - row_top * cell_size_y_map
        )
        
        feedback.setProgressText(f"Reading raster data ({read_width} × {read_height} cells covering polygons)...")
        block = dp.block(1, read_extent_rect, read_width, read_height)
        
        if block is None or block.isEmpty():
            raise QgsProcessingException("Failed to read raster data")

        # Extract raster values (will filter by polygons during processing)
        raster_data = {}  # (col, row) -> value
        cell_values = {}  # (col, row) -> (x, y, value)
        
        for row_idx in range(read_height):
            if feedback.isCanceled():
                return {self.OUTPUT: dest_id}
            for col_idx in range(read_width):
                try:
                    val = block.value(row_idx, col_idx)
                    if val is not None:
                        try:
                            zv = float(val)
                            if nodata is None or abs(zv - nodata) > 1e-9:
                                col = read_col_start + col_idx
                                row = read_row_start + row_idx
                                if 0 <= col < width and 0 <= row < height:
                                    x, y = cell_centroid_to_map(col, row)
                                    raster_data[(col, row)] = zv
                                    cell_values[(col, row)] = (x, y, zv)
                        except (ValueError, TypeError):
                            pass
                except Exception:
                    pass

        if not raster_data:
            raise QgsProcessingException("No valid raster cells found in processing area")

        # Process each polygon separately
        all_points = []
        
        feedback.setProgressText(f"Finding extrema within {len(polygon_features)} polygon(s)...")
        
        for poly_idx, (poly_geom, poly_feat) in enumerate(polygon_features):
            if feedback.isCanceled():
                break
            
            # Find all raster cells within this polygon
            polygon_cells = []
            for (col, row), (x, y, val) in cell_values.items():
                cell_point = QgsPointXY(x, y)
                if poly_geom.contains(QgsGeometry.fromPointXY(cell_point)):
                    polygon_cells.append((x, y, val))
            
            if not polygon_cells:
                continue
            
            # Find extrema within this polygon
            values = [c[2] for c in polygon_cells]
            
            if find_maxima:
                max_val = max(values)
                max_cells = [c for c in polygon_cells if abs(c[2] - max_val) < 1e-9]
                if max_cells:
                    # Use centroid of all cells with max value
                    cx = sum(c[0] for c in max_cells) / len(max_cells)
                    cy = sum(c[1] for c in max_cells) / len(max_cells)
                    all_points.append({
                        'x': cx,
                        'y': cy,
                        'type': 'maximum',
                        'value': max_val,
                        'polygon_id': poly_feat.id()
                    })
            
            if find_minima:
                min_val = min(values)
                min_cells = [c for c in polygon_cells if abs(c[2] - min_val) < 1e-9]
                if min_cells:
                    # Use centroid of all cells with min value
                    cx = sum(c[0] for c in min_cells) / len(min_cells)
                    cy = sum(c[1] for c in min_cells) / len(min_cells)
                    all_points.append({
                        'x': cx,
                        'y': cy,
                        'type': 'minimum',
                        'value': min_val,
                        'polygon_id': poly_feat.id()
                    })
            
            feedback.setProgress(int(50 * (poly_idx + 1) / len(polygon_features)))

        # Write output points
        feedback.setProgressText(f"Writing {len(all_points)} output points...")
        for idx, pt in enumerate(all_points):
            if feedback.isCanceled():
                break
            feat = QgsFeature(fields)
            feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(pt['x'], pt['y'])))
            feat.setAttributes([pt['type'], pt['value'], pt['polygon_id']])
            sink.addFeature(feat, QgsFeatureSink.FastInsert)
            
            if all_points:
                feedback.setProgress(int(50 + 50 * (idx + 1) / len(all_points)))

        # Style output layer
        layer = context.temporaryLayerStore().mapLayer(dest_id) or context.getMapLayer(dest_id)
        if layer:
            root = QgsRuleBasedRenderer.Rule(None)
            
            if find_maxima:
                rule_max = QgsRuleBasedRenderer.Rule(
                    QgsMarkerSymbol.createSimple({
                        "color": "#cc3333",
                        "outline_color": "#cc3333",
                        "size": "2",
                        "outline_width": "0",
                        "name": "circle"
                    })
                )
                rule_max.setFilterExpression('"extrema_type" = \'maximum\'')
                rule_max.setLabel("Maximum")
                root.appendChild(rule_max)
            
            if find_minima:
                rule_min = QgsRuleBasedRenderer.Rule(
                    QgsMarkerSymbol.createSimple({
                        "color": "#2ca02c",
                        "outline_color": "#2ca02c",
                        "size": "2",
                        "outline_width": "0",
                        "name": "circle"
                    })
                )
                rule_min.setFilterExpression('"extrema_type" = \'minimum\'')
                rule_min.setLabel("Minimum")
                root.appendChild(rule_min)
            
            layer.setRenderer(QgsRuleBasedRenderer(root))
            layer.triggerRepaint()

        return {self.OUTPUT: dest_id}

