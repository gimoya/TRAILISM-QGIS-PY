"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterString,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFile,
    QgsProcessingException,
    QgsPoint,
    QgsPointXY,
    QgsGeometry,
    QgsWkbTypes,
    QgsVectorLayer,
    QgsProject,
)
import math
import json
import struct
import base64
from pathlib import Path

# Constants for curve detection and simple shapes
MIN_SHAPE_HEIGHT = 0.25  # Minimum height value for simple shapes
DENSIFY_DISTANCE = 0.25  # Fixed densification distance (meters) for curve detection and geometry processing
CURVE_WINDOW_DISTANCE_DEFAULT = 8 # Default curve detection window distance (meters)
CURVE_THRESHOLD_ANGLE = 10  # Single threshold for curve detection (degrees) - average of previous start/end
CURVE_WINDOW_FORWARD_RATIO = 0.4  # 40% of window is forward (for late start)
CURVE_WINDOW_BACKWARD_RATIO = 0.6  # 60% of window is backward (for late end)
MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION = 10  # Minimum curve length (meters) for full banking height. Shorter curves get reduced height.
CURVE_LENGTH_REDUCTION_RATIO = 0.05  # Height reduction ratio per meter for curves below MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION (0.1 = 10% per meter)

class AddCurveBankingToLineAlgorithm(QgsProcessingAlgorithm):
    """
    Detects curves in a linestring by analyzing bisector angles over a distance window.
    For detected curves, widens the ribbon and raises only the outer edge of the curve.
    """

    INPUT = 'INPUT'
    DEM = 'DEM'
    SELECTED_FEATURES_ONLY = 'SELECTED_FEATURES_ONLY'
    Z_RAISE_PARAM = 'Z_RAISE'
    CURVE_WINDOW_DISTANCE = 'CURVE_WINDOW_DISTANCE'
    OUTER_EDGE_RAISE = 'OUTER_EDGE_RAISE'
    TRANSVERSAL_CONCAVITY = 'TRANSVERSAL_CONCAVITY'
    RIBBON_WIDTH = 'RIBBON_WIDTH'
    EXPORT_FORMAT = 'EXPORT_FORMAT'
    TEXTURE_FILE = 'TEXTURE_FILE'
    ONLY_RAISE_TRAIL = 'ONLY_RAISE_TRAIL'
    SHAPE_DISTANCES = 'SHAPE_DISTANCES'
    SHAPE_HEIGHTS = 'SHAPE_HEIGHTS'
    SHAPE_TABLE_LENGTHS = 'SHAPE_TABLE_LENGTHS'
    SIGMOID_LENGTHS = 'SIGMOID_LENGTHS'
    SIGMOID_STEEPNESS_K_PARAM = 'SIGMOID_STEEPNESS_K'
    SIMPLE_SHAPES = 'SIMPLE_SHAPES'
    DENSIFY_DISTANCE = 'DENSIFY_DISTANCE'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return AddCurveBankingToLineAlgorithm()

    def name(self):
        return 'add_trail_shapes'

    def displayName(self):
        return 'TRAILISM: Add Trail Shapes'
    
    def group(self):
        return 'TRAILISM'
    
    def groupId(self):
        return 'trailism'

    def shortHelpString(self):
        return (
            'TRAILISM ALGORITHM:\n\n'
            'IMPORTANT: All numeric values use period (.) as decimal separator for consistency.\n'
            'Comma (,) is used only to separate list items.\n'
            'Note: The GUI may display default values with comma (,) due to system locale, but always enter values using period (.) as decimal separator.\n\n'
            'Detects curves in a linestring by analyzing bisector angles over a distance window.\n'
            'For detected curves, raises only the outer edge of the curve with a smooth parabolic arc.\n'
            'Optionally adds smooth sigmoid "shape" bumps at specified distances.\n'
            'Shapes and curve banking can be combined for complex trail profiles.\n\n'
            'Parameters:\n'
            '- DEM layer: Required raster layer for sampling base Z values\n'
            '- Trail raise height: Height to raise the trail above terrain\n'
            '  Default: 1.0 meters, minimum: 0.0\n'
            '  Use period (.) for decimal separator\n'
            '- Outer edge raise: Height to raise the outer edge of curves (meters)\n'
            '  Default: 0.3 meters, minimum: 0.0\n'
            '  Use period (.) for decimal separator\n'
            '- Ribbon width: Base width of the exported ribbon mesh in meters\n'
            '  Default: 2.0 meters, minimum: 0.5\n'
            '  Use period (.) for decimal separator\n'
            '- Only raise trail: If checked, skip shapes and curve banking, just raise trail\n'
            '- Simple shapes mode: If checked, only distances and heights are required.\n'
            '  Table and sigmoid lengths are auto-calculated: table = h*3 - 3 * minimum height, sigmoid = h*4\n'
            '  Minimum height: 0.25 meters\n'
            '- Shape distances: Comma-separated distances from start (e.g., "10,30,end")\n'
            '  Use period (.) for decimal values, "end" for line end, "0" for start\n'
            '- Shape heights: Comma-separated heights for each shape (must match distances count)\n'
            '  Minimum: 0.25 meters (enforced in simple mode)\n'
            '- Shape table lengths: Comma-separated flat portion lengths (must match distances count)\n'
            '  Required only if Simple shapes mode is unchecked\n'
            '- Sigmoid transition lengths: Comma-separated lengths for sigmoid curves (must match distances count)\n'
            '  Required only if Simple shapes mode is unchecked\n'
            '- Sigmoid steepness factor K: Controls transition steepness (1.0 to 30.0, default: 4.5)\n'
            '- Densify distance: Maximum distance between vertices (default: 0.25m, minimum: 0.05m)\n'
            '- Export format: GLTF or GLB format for 3D visualization\n\n'
            'Output:\n'
            '- Exports ribbon meshes to GLTF/GLB format.\n'
            '- Original layer is not modified.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT,
                'Input line layer',
                [QgsProcessing.TypeVectorLine]
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.ONLY_RAISE_TRAIL,
                'Only raise trail, skip shapes and curve banking',
                defaultValue=False
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.SELECTED_FEATURES_ONLY,
                'Selected features only',
                defaultValue=True
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.DEM,
                'DEM layer',
                optional=False
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.Z_RAISE_PARAM,
                'Trail raise height (meters). Default: 0.5, minimum: 0.1. Use period (.) for decimal separator',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
                minValue=0.1
            )
        )

        # Curve detection parameters are now internal constants (not exposed in GUI)
        # CURVE_WINDOW_DISTANCE_DEFAULT = 15.0
        # CURVE_THRESHOLD_ANGLE_DEFAULT = 10.0

        self.addParameter(
            QgsProcessingParameterNumber(
                self.OUTER_EDGE_RAISE,
                'Outer edge raise (meters). Default: 1, minimum: 0.0. Use period (.) for decimal separator',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1,
                minValue=0.0,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.TRANSVERSAL_CONCAVITY,
                'Transversal concavity depth (meters). Default: 0.3, positive = concave. Range: 0.0 to 1.0. Use period (.) for decimal separator',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.3,
                minValue=0.0,
                maxValue=1.0,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.RIBBON_WIDTH,
                'Ribbon width (meters). Default: 2.0, minimum: 0.5. Use period (.) for decimal separator',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=2,
                minValue=0.5,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterFile(
                self.TEXTURE_FILE,
                'Texture image file (optional)',
                behavior=QgsProcessingParameterFile.File,
                extension='png,jpg,jpeg',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.SIMPLE_SHAPES,
                'Simple shapes mode (auto-calculate table and sigmoid lengths from heights)',
                defaultValue=False
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_DISTANCES,
                'Shape distances from start (comma-separated, meters). Use period (.) for decimals. Use "end" for line end',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_HEIGHTS,
                'Shape heights (comma-separated, meters). Use period (.) for decimals. Minimum: 0.25',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_TABLE_LENGTHS,
                'Shape table lengths (comma-separated, meters) - flat portion between sigmoid start and end. Use period (.) for decimals',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SIGMOID_LENGTHS,
                'Sigmoid transition lengths (comma-separated, meters) - length of sigmoid curve at start and end of each shape. Must match number of shape distances. Use period (.) for decimals',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.SIGMOID_STEEPNESS_K_PARAM,
                'Sigmoid steepness factor K (higher = steeper transition, lower = flatter transition). Range: 1.0 to 30.0. Default: 4.5. Use period (.) for decimal separator',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=4.5,
                minValue=1,
                maxValue=30,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.EXPORT_FORMAT,
                'Export format',
                options=['GLTF', 'GLB'],
                defaultValue=1  # GLB by default
            )
        )


    def processAlgorithm(self, parameters, context, feedback):
        # Get input layer directly (must be a vector layer loaded in project)
        input_layer = self.parameterAsVectorLayer(parameters, self.INPUT, context)
        if input_layer is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))
        
        # Verify it's a vector layer
        if not isinstance(input_layer, QgsVectorLayer):
            raise QgsProcessingException(
                f"Layer '{input_layer.name()}' is not a vector layer."
            )
        
        # Check "only raise trail" FIRST - if enabled, skip all other parameters
        only_raise_trail = self.parameterAsBool(parameters, self.ONLY_RAISE_TRAIL, context)
        
        dem_layer = self.parameterAsRasterLayer(parameters, self.DEM, context)
        if dem_layer is None:
            raise QgsProcessingException("DEM layer is required")
        
        # Read basic parameters (always needed)
        z_raise = self.parameterAsDouble(parameters, self.Z_RAISE_PARAM, context)
        z_raise = float(z_raise) if z_raise is not None else 0.5
        if z_raise < 0.1:
            z_raise = 0.1
        
        # Use fixed densification distance from constant
        densify_distance = DENSIFY_DISTANCE
        ribbon_width = self.parameterAsDouble(parameters, self.RIBBON_WIDTH, context)
        ribbon_width = float(ribbon_width) if ribbon_width is not None else 3.0
        export_format_idx = self.parameterAsEnum(parameters, self.EXPORT_FORMAT, context)
        export_format = 'GLB' if export_format_idx == 1 else 'GLTF'
        texture_file = self.parameterAsString(parameters, self.TEXTURE_FILE, context) or None
        
        # Skip all curve banking and shape parameters if "only raise trail" is enabled
        if only_raise_trail:
            curve_window_distance = 0.0
            curve_threshold_angle = 0.0
            outer_edge_raise = 0.0
            transversal_concavity = 0.0
            simple_shapes = False
            shape_distances_str = None
            shape_heights_str = None
            shape_table_lengths_str = None
            sigmoid_lengths_str = None
            sigmoid_steepness_k = None
        else:
            # Read curve banking parameters
            # Use internal constants for curve detection (not exposed as parameters)
            curve_window_distance = CURVE_WINDOW_DISTANCE_DEFAULT
            # Use single threshold for curve detection
            curve_threshold_angle = CURVE_THRESHOLD_ANGLE
            outer_edge_raise = self.parameterAsDouble(parameters, self.OUTER_EDGE_RAISE, context)
            outer_edge_raise = float(outer_edge_raise) if outer_edge_raise is not None else 1.0
            transversal_concavity = self.parameterAsDouble(parameters, self.TRANSVERSAL_CONCAVITY, context)
            transversal_concavity = float(transversal_concavity) if transversal_concavity is not None else 0.0
            
            # Read shape parameters
            simple_shapes = self.parameterAsBool(parameters, self.SIMPLE_SHAPES, context)
            shape_distances_str = self.parameterAsString(parameters, self.SHAPE_DISTANCES, context) or ''
            shape_heights_str = self.parameterAsString(parameters, self.SHAPE_HEIGHTS, context) or ''
            shape_table_lengths_str = self.parameterAsString(parameters, self.SHAPE_TABLE_LENGTHS, context) or ''
            sigmoid_lengths_str = self.parameterAsString(parameters, self.SIGMOID_LENGTHS, context) or ''
            sigmoid_steepness_k = self.parameterAsDouble(parameters, self.SIGMOID_STEEPNESS_K_PARAM, context)
            sigmoid_steepness_k = float(sigmoid_steepness_k) if sigmoid_steepness_k is not None else 4.5
        
        if feedback:
            feedback.pushInfo(f"Trail raise height: {z_raise} meters")
            feedback.pushInfo(f"Only raise trail: {only_raise_trail}")
            if not only_raise_trail:
                feedback.pushInfo(f"Curve detection window: {curve_window_distance} meters (40% forward, 60% backward)")
                feedback.pushInfo(f"Curve threshold angle: {CURVE_THRESHOLD_ANGLE}°")
                feedback.pushInfo(f"Outer edge raise: {outer_edge_raise} meters")
                feedback.pushInfo(f"Transversal concavity: {transversal_concavity} meters")
        
        # Get project home directory for ribbon export
        project = QgsProject.instance()
        project_home = project.homePath()
        if not project_home or not Path(project_home).exists():
            project_home = Path(__file__).parent.absolute()
        else:
            project_home = Path(project_home)
        
        # Get DEM provider
        dem_provider = dem_layer.dataProvider()
        if dem_provider is None:
            raise QgsProcessingException("Invalid DEM layer")
        
        # Check if only selected features should be processed
        selected_only = self.parameterAsBool(parameters, self.SELECTED_FEATURES_ONLY, context)
        
        # Get features to process
        if selected_only:
            selected_ids = input_layer.selectedFeatureIds()
            if not selected_ids:
                raise QgsProcessingException("No features are selected. Please select features or uncheck 'Selected features only'.")
            features = [f for f in input_layer.getFeatures(selected_ids)]
            feedback.pushInfo(f"Processing {len(features)} selected features for ribbon export...")
        else:
            features = list(input_layer.getFeatures())
            feedback.pushInfo(f"Processing {len(features)} features for ribbon export...")
        
        if not features:
            raise QgsProcessingException("No features to process.")
        
        # Calculate progress increment once
        total = 100.0 / len(features) if features else 0
        
        # Track which features are exported
        exported_features = []
        
        try:
            for current, feature in enumerate(features):
                if feedback.isCanceled():
                    break
                
                geom = feature.geometry()
                if not geom or geom.isEmpty():
                    feedback.pushWarning(f"Feature {feature.id()} has empty or invalid geometry, skipping.")
                    continue
                
                if geom.type() != QgsWkbTypes.LineGeometry:
                    feedback.pushWarning(f"Feature {feature.id()} is not a line geometry, skipping.")
                    continue
                
                # Note: Z values from input geometry are ignored - only XY coordinates are used
                # All Z values come from DEM sampling
                
                try:
                    if only_raise_trail:
                        # Simple path: just raise trail, no shapes or curve banking
                        result = self.apply_raise_only(
                            geom,
                            dem_provider,
                            z_raise,
                            feedback
                        )
                        
                        if isinstance(result, tuple):
                            new_geom, base_z_values, smoothed_z_values = result
                        else:
                            new_geom = result
                            base_z_values = None
                            smoothed_z_values = None
                        
                        if new_geom is None or new_geom.isEmpty():
                            feedback.pushWarning(f"Feature {feature.id()}: Failed to create output geometry, skipping.")
                            continue
                        
                        # Export ribbon mesh (simple, no curve banking)
                        if base_z_values is not None and smoothed_z_values is not None:
                            try:
                                self.export_ribbon_mesh_simple(
                                    feature,
                                    new_geom,
                                    base_z_values,
                                    smoothed_z_values,
                                    ribbon_width,
                                    export_format,
                                    project_home,
                                    texture_file,
                                    input_layer,
                                    dem_layer,
                                    z_raise,
                                    feedback
                                )
                                exported_features.append(feature.id())
                            except Exception as e:
                                feedback.pushWarning(f"Feature {feature.id()}: Failed to export ribbon mesh - {str(e)}")
                    else:
                        # Parse shape parameters if provided
                        shape_distances_raw = []
                        shape_heights = []
                        shape_table_lengths = []
                        sigmoid_lengths = []
                        first_is_zero = False
                        last_is_end = False
                        use_shapes = False
                        
                        def parse_float_list(value_str):
                            items = []
                            for item in value_str.split(','):
                                item = item.strip()
                                if item:
                                    items.append(item)
                            return items
                        
                        # Check if shapes should be used
                        if simple_shapes:
                            # Simple mode: only need distances and heights
                            if shape_distances_str and shape_heights_str:
                                try:
                                    shape_distances_raw = parse_float_list(shape_distances_str)
                                    heights_parsed = parse_float_list(shape_heights_str)
                                    
                                    if len(shape_distances_raw) == len(heights_parsed) and len(shape_distances_raw) > 0:
                                        # Validate and parse heights (minimum MIN_SHAPE_HEIGHT)
                                        shape_heights = []
                                        for h_str in heights_parsed:
                                            h = float(h_str)
                                            if h < MIN_SHAPE_HEIGHT:
                                                if feedback:
                                                    feedback.pushWarning(f"Feature {feature.id()}: Shape height {h} is below minimum {MIN_SHAPE_HEIGHT}, using {MIN_SHAPE_HEIGHT}")
                                                h = MIN_SHAPE_HEIGHT
                                            shape_heights.append(h)
                                        
                                        # Calculate table lengths and sigmoid lengths from heights
                                        # table_length = h * 3 - min_height * 3 
                                        # sigmoid_length = h * 4
                                        shape_table_lengths = [h * 3.0 - MIN_SHAPE_HEIGHT * 3.0 for h in shape_heights]
                                        sigmoid_lengths = [h * 4.0 for h in shape_heights]
                                        
                                        use_shapes = True
                                        
                                        # Calculate line length for 'end' keyword
                                        line_length = geom.length()
                                        if line_length > 0:
                                            # Parse distances, replacing 'end' with actual line length
                                            for idx, dist_str in enumerate(shape_distances_raw):
                                                dist_str_lower = dist_str.lower().strip()
                                                if dist_str_lower == 'end':
                                                    if idx == len(shape_distances_raw) - 1:
                                                        last_is_end = True
                                                    shape_distances_raw[idx] = str(line_length)
                                                elif dist_str_lower == '0' or dist_str == '0':
                                                    if idx == 0:
                                                        first_is_zero = True
                                        
                                        if feedback:
                                            feedback.pushInfo(f"Feature {feature.id()}: Simple shapes mode - calculated {len(shape_heights)} shapes from heights")
                                            for i, (h, t, s) in enumerate(zip(shape_heights, shape_table_lengths, sigmoid_lengths)):
                                                feedback.pushInfo(f"  Shape {i+1}: height={h:.2f}m, table_length={t:.2f}m, sigmoid_length={s:.2f}m")
                                except Exception as e:
                                    if feedback:
                                        feedback.pushWarning(f"Feature {feature.id()}: Error parsing simple shape parameters - {str(e)}, continuing without shapes")
                                    use_shapes = False
                        else:
                            # Advanced mode: require all parameters
                            if shape_distances_str and shape_heights_str and shape_table_lengths_str and sigmoid_lengths_str:
                                try:
                                    shape_distances_raw = parse_float_list(shape_distances_str)
                                    shape_heights = [float(h) for h in parse_float_list(shape_heights_str)]
                                    shape_table_lengths = [float(t) for t in parse_float_list(shape_table_lengths_str)]
                                    sigmoid_lengths = [float(s) for s in parse_float_list(sigmoid_lengths_str)]
                                    
                                    if (len(shape_distances_raw) == len(shape_heights) == len(shape_table_lengths) == len(sigmoid_lengths) and 
                                        len(shape_distances_raw) > 0):
                                        use_shapes = True
                                        
                                        # Calculate line length for 'end' keyword
                                        line_length = geom.length()
                                        if line_length > 0:
                                            # Parse distances, replacing 'end' with actual line length
                                            for idx, dist_str in enumerate(shape_distances_raw):
                                                dist_str_lower = dist_str.lower().strip()
                                                if dist_str_lower == 'end':
                                                    if idx == len(shape_distances_raw) - 1:
                                                        last_is_end = True
                                                    shape_distances_raw[idx] = str(line_length)
                                                elif dist_str_lower == '0' or dist_str == '0':
                                                    if idx == 0:
                                                        first_is_zero = True
                                except Exception as e:
                                    if feedback:
                                        feedback.pushWarning(f"Feature {feature.id()}: Error parsing shape parameters - {str(e)}, continuing without shapes")
                                    use_shapes = False
                        
                        # Apply curve banking (and shapes if provided)
                        # Since we only use length-based reduction, threshold is only for validation
                        curve_threshold_angle_rad = math.radians(CURVE_THRESHOLD_ANGLE) if not only_raise_trail else 0.0
                        
                    result = self.apply_curve_banking(
                        geom,
                        dem_provider,
                        z_raise,
                            curve_window_distance if not only_raise_trail else 0.0,
                        curve_threshold_angle_rad,
                            outer_edge_raise if not only_raise_trail else 0.0,
                            use_shapes,
                            shape_distances_raw if use_shapes else [],
                            shape_heights if use_shapes else [],
                            shape_table_lengths if use_shapes else [],
                            sigmoid_lengths if use_shapes else [],
                            sigmoid_steepness_k if use_shapes else 4.5,
                            first_is_zero,
                            last_is_end,
                        feedback
                    )
                    
                    if isinstance(result, tuple) and len(result) >= 8:
                        original_geom, base_z_values, smoothed_z_values, curve_flags, curve_directions, banking_heights, distances, shape_z_offsets = result
                    elif isinstance(result, tuple) and len(result) >= 7:
                        # Fallback for old return format
                        original_geom, base_z_values, smoothed_z_values, curve_flags, curve_directions, banking_heights, distances = result
                        shape_z_offsets = [0.0] * len(distances) if distances else []
                    else:
                        feedback.pushWarning(f"Feature {feature.id()}: Failed to process, skipping.")
                        continue
                    
                    if original_geom is None or original_geom.isEmpty():
                        feedback.pushWarning(f"Feature {feature.id()}: Failed to create output geometry, skipping.")
                        continue
                    
                    # Export ribbon mesh
                    if base_z_values is not None and smoothed_z_values is not None:
                        try:
                            self.export_ribbon_mesh(
                                feature,
                                original_geom,
                                base_z_values,
                                smoothed_z_values,
                                curve_flags,
                                curve_directions,
                                banking_heights,
                                distances,
                                shape_z_offsets,
                                ribbon_width,
                                transversal_concavity,
                                export_format,
                                project_home,
                                texture_file,
                                input_layer,
                                dem_layer,
                                z_raise,
                                feedback
                            )
                            exported_features.append(feature.id())
                        except Exception as e:
                            feedback.pushWarning(f"Feature {feature.id()}: Failed to export ribbon mesh - {str(e)}")
                except Exception as e:
                    feedback.pushWarning(f"Feature {feature.id()}: Error processing - {str(e)}, skipping.")
                    continue
                
                feedback.setProgress(int(current * total))
            
            # Report export results
            if exported_features:
                feedback.pushInfo(f"Processing complete. Ribbon meshes exported for {len(exported_features)} feature(s): {', '.join(map(str, exported_features))}")
            else:
                feedback.pushInfo("Processing complete. No ribbon meshes were exported.")
        except Exception as e:
            raise QgsProcessingException(f"Error processing layer: {str(e)}")
        
        return {}
    
    @staticmethod
    def normalized_sigmoid(t, k):
        """
        Calculate normalized sigmoid value that goes exactly from 0 to 1.
        This ensures the transition completes within the specified length regardless of steepness.
        
        Args:
            t: Parameter from 0 to 1
            k: Steepness factor (higher = steeper, lower = flatter)
        
        Returns:
            Normalized sigmoid value from 0 to 1
        """
        # Calculate raw sigmoid
        raw_sigmoid = 1.0 / (1.0 + math.exp(-k * (2.0 * t - 1.0)))
        
        # Calculate sigmoid values at boundaries
        sigmoid_at_0 = 1.0 / (1.0 + math.exp(k))
        sigmoid_at_1 = 1.0 / (1.0 + math.exp(-k))
        sigmoid_range = sigmoid_at_1 - sigmoid_at_0
        
        # Normalize to exactly 0-1 range
        if sigmoid_range > 1e-9:
            normalized = (raw_sigmoid - sigmoid_at_0) / sigmoid_range
        else:
            normalized = raw_sigmoid
        
        return max(0.0, min(1.0, normalized))  # Clamp to [0, 1]
    
    def apply_raise_only(self, geometry, dem_provider, z_raise=1.0, feedback=None):
        """
        Apply only Z raise to trail, no shapes or curve banking.
        Drapes from DEM, smooths Z values, then adds raise offset.
        """
        # Validate geometry
        if geometry is None or geometry.isEmpty():
            return (geometry, None, None)
        
        # Densify geometry to exact intervals if requested
        # Always densify using fixed distance
        geometry = self._densify_at_exact_intervals(geometry, DENSIFY_DISTANCE, feedback)
        
        # Get vertices
        try:
            vertices = list(geometry.vertices())
        except Exception:
            return (geometry, None, None)
        
        if len(vertices) < 2:
            return (geometry, None, None)
        
        # Calculate cumulative distances along the line
        distances = [0.0]
        total_length = 0.0
        
        for i in range(1, len(vertices)):
            prev_vertex = vertices[i-1]
            curr_vertex = vertices[i]
            
            dx = curr_vertex.x() - prev_vertex.x()
            dy = curr_vertex.y() - prev_vertex.y()
            segment_length = math.sqrt(dx*dx + dy*dy)
            
            total_length += segment_length
            distances.append(total_length)
        
        # Step 1: Drape Z values from DEM at all vertices
        base_z_values = []
        for i, vertex in enumerate(vertices):
            x = vertex.x()
            y = vertex.y()
            val, ok = dem_provider.sample(QgsPointXY(x, y), 1)
            if ok and val is not None:
                base_z_values.append(float(val))
            else:
                base_z_values.append(0.0)
        
        # Step 2: Smooth the base Z values using 10 meter distance-based windows
        smoothing_distance = 10.0
        half_window = smoothing_distance / 2.0
        smoothed_z_values = []
        
        num_vertices = len(base_z_values)
        window_start_idx = 0
        
        for i in range(num_vertices):
            current_distance = distances[i]
            window_start_dist = current_distance - half_window
            window_end_dist = current_distance + half_window
            
            while window_start_idx < num_vertices and distances[window_start_idx] < window_start_dist:
                window_start_idx += 1
            
            window_z_values = []
            j = window_start_idx
            while j < num_vertices and distances[j] <= window_end_dist:
                window_z_values.append(base_z_values[j])
                j += 1
            
            if window_z_values:
                smoothed_z = sum(window_z_values) / len(window_z_values)
            else:
                smoothed_z = base_z_values[i]
            
            smoothed_z_values.append(smoothed_z)
        
        # Step 3: Create new points with only raise added
        if z_raise is None:
            z_raise = 0.0
        z_raise = float(z_raise)
        
        new_points = []
        for i, vertex in enumerate(vertices):
            x = vertex.x()
            y = vertex.y()
            base_z = smoothed_z_values[i]
            z = base_z + z_raise
            new_points.append(QgsPoint(x, y, z))
        
        # Recreate geometry
        new_geometry = QgsGeometry.fromPolyline(new_points)
        
        return (new_geometry, base_z_values, smoothed_z_values)
    
    def create_parabolic_arc_envelope(self, curve_flags, distances, transition_points,
                                      outer_edge_raise, curve_directions=None, feedback=None):
        """
        Create smooth parabolic arc envelope for banking.
        Rises smoothly, peaks in the middle, falls smoothly - like a regular arc.
        Uses full curve segment for envelope to ensure smooth transitions between closely spaced curves.
        Banking height is length-based only (no angle-based scaling).
        Blends transitions between curves that are within 5m of each other.
        Adds transition zones where banking fades to zero when curve direction changes (left to right or vice versa).
        Short curves (below MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION) get reduced banking height based on meters below threshold.
        """
        banking_heights = [0.0] * len(curve_flags)
        if curve_directions is None:
            curve_directions = [0.0] * len(curve_flags)
        if transition_points is None:
            transition_points = set()
        
        # First pass: identify all curve segments with their directions and distances
        curve_segments_info = []
        in_curve = False
        curve_start_idx = None
        
        for i in range(len(curve_flags)):
            if curve_flags[i] and not in_curve:
                # Start of new curve segment
                in_curve = True
                curve_start_idx = i
            elif not curve_flags[i] and in_curve:
                # End of curve segment
                in_curve = False
                curve_end_idx = i - 1
                # Get direction for this segment
                mid_idx = (curve_start_idx + curve_end_idx) // 2
                seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
                curve_segments_info.append({
                    'start_idx': curve_start_idx,
                    'end_idx': curve_end_idx,
                    'start_dist': distances[curve_start_idx],
                    'end_dist': distances[curve_end_idx],
                    'direction': seg_dir
                })
        
        if in_curve:
            curve_end_idx = len(curve_flags) - 1
            mid_idx = (curve_start_idx + curve_end_idx) // 2
            seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
            curve_segments_info.append({
                'start_idx': curve_start_idx,
                'end_idx': curve_end_idx,
                'start_dist': distances[curve_start_idx],
                'end_dist': distances[curve_end_idx],
                'direction': seg_dir
            })
        
        # Create lookup: which segment does each index belong to?
        idx_to_seg = {}
        for seg_idx, seg in enumerate(curve_segments_info):
            for j in range(seg['start_idx'], seg['end_idx'] + 1):
                idx_to_seg[j] = seg_idx
        
        # Second pass: create envelope, checking for opposite-direction neighbors
        in_curve = False
        curve_start_idx = None
        
        for i in range(len(curve_flags)):
            if curve_flags[i] and not in_curve:
                # Start of new curve segment
                in_curve = True
                curve_start_idx = i
            elif not curve_flags[i] and in_curve:
                # End of curve segment
                in_curve = False
                curve_end_idx = i - 1
                
                # Clamp curve boundaries to transition points FIRST to prevent banking from crossing boundaries
                # Find nearest transition points to clamp boundaries
                clamped_start_idx = curve_start_idx
                clamped_end_idx = curve_end_idx
                
                # Find last transition point before or at start
                for tp in sorted(transition_points, reverse=True):
                    if tp < curve_start_idx:
                        clamped_start_idx = tp + 1  # Start after transition point
                        break
                    elif tp == curve_start_idx:
                        # Transition point exactly at start - start after it
                        clamped_start_idx = tp + 1
                        break
                
                # Find first transition point at or after end, or within the curve
                for tp in sorted(transition_points):
                    if tp > curve_end_idx:
                        clamped_end_idx = tp - 1  # End before transition point
                        break
                    elif tp >= curve_start_idx and tp <= curve_end_idx:
                        # Transition point inside curve - clamp end to before it
                        clamped_end_idx = tp - 1
                        break
                
                # Ensure valid range
                clamped_start_idx = max(0, clamped_start_idx)
                clamped_end_idx = min(len(curve_flags) - 1, clamped_end_idx)
                
                # Only apply banking if clamped range is valid
                if clamped_start_idx > clamped_end_idx:
                    continue
                
                # Use clamped boundaries for all calculations
                curve_start_dist = distances[clamped_start_idx]
                curve_end_dist = distances[clamped_end_idx]
                curve_length = curve_end_dist - curve_start_dist
                
                if feedback:
                    dir_str = "LEFT" if seg_dir > 0 else "RIGHT"
                    feedback.pushInfo(f"Banking envelope {dir_str}: vertices {clamped_start_idx}-{clamped_end_idx}, distances {curve_start_dist:.2f}-{curve_end_dist:.2f}m, length {curve_length:.2f}m")
                
                if curve_length <= 0:
                    continue
                
                # Get segment direction BEFORE calculating peak height
                current_seg_idx = idx_to_seg.get(curve_start_idx, -1)
                seg_dir = curve_segments_info[current_seg_idx]['direction'] if current_seg_idx >= 0 else 0.0
                
                # Calculate peak height for clamped segment
                peak_height = outer_edge_raise
                if curve_length < MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION:
                    meters_below_threshold = MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION - curve_length
                    reduction_factor = max(0.0, 1.0 - meters_below_threshold * CURVE_LENGTH_REDUCTION_RATIO)
                    peak_height *= reduction_factor
                
                # Determine fade zones at boundaries (based on original segment boundaries, not clamped)
                # This applies to ALL curves, not just short ones
                fade_zone = 5.0  # meters
                transition_distance = 5.0  # meters - minimum gap for fade
                fade_start = 0.0
                fade_end = 0.0
                
                # Check previous segment
                if current_seg_idx > 0:
                    prev_seg = curve_segments_info[current_seg_idx - 1]
                    gap = curve_start_dist - prev_seg['end_dist']
                    # Only fade if there's a straight section (gap >= transition_distance)
                    # For back-to-back curves (gap < transition_distance), keep full banking
                    if gap >= transition_distance:
                        prev_dir = prev_seg['direction']
                        if (prev_dir > 0 and seg_dir < 0) or (prev_dir < 0 and seg_dir > 0):
                            # Opposite direction - fade in at start
                            fade_start = min(fade_zone, curve_length / 2.0)
                
                # Check next segment
                if current_seg_idx < len(curve_segments_info) - 1:
                    next_seg = curve_segments_info[current_seg_idx + 1]
                    gap = next_seg['start_dist'] - curve_end_dist
                    # Only fade if there's a straight section (gap >= transition_distance)
                    # For back-to-back curves (gap < transition_distance), keep full banking
                    if gap >= transition_distance:
                        next_dir = next_seg['direction']
                        if (seg_dir > 0 and next_dir < 0) or (seg_dir < 0 and next_dir > 0):
                            # Opposite direction - fade out at end
                            fade_end = min(fade_zone, curve_length / 2.0)
                    
                # Calculate banking for segment vertices (using clamped range)
                # Envelope is recalculated over the clamped range to ensure smooth transitions
                # This applies to ALL curves, not just short ones
                for j in range(clamped_start_idx, clamped_end_idx + 1):
                    # Position within clamped curve: 0.0 (clamped start) to 1.0 (clamped end)
                    # This ensures envelope starts at 0 and ends at 0 within the clamped range
                    # Clamp t to [0, 1] to handle floating point precision issues at boundaries
                    t = (distances[j] - curve_start_dist) / curve_length
                    t = max(0.0, min(1.0, t))  # Clamp to [0, 1]
                    
                    # Parabolic arc: y = 4 * x * (1 - x)
                    # This gives a smooth curve that is 0 at t=0 and t=1, peaking at t=0.5
                    envelope = 4.0 * t * (1.0 - t)
                    
                    # Apply fade at boundaries if opposite directions (only if boundary wasn't clamped)
                    # Clamped boundaries already fade to 0 naturally via the envelope
                    was_clamped_start = clamped_start_idx != curve_start_idx
                    was_clamped_end = clamped_end_idx != curve_end_idx
                    
                    if not was_clamped_start and fade_start > 0 and t < fade_start / curve_length:
                        fade_t = t * curve_length / fade_start
                        fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)  # smoothstep
                        envelope *= fade_factor
                    
                    if not was_clamped_end and fade_end > 0 and t > 1.0 - (fade_end / curve_length):
                        fade_t = (1.0 - t) * curve_length / fade_end
                        fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)  # smoothstep
                        envelope *= fade_factor
                    
                    # Blend with existing banking if curves overlap
                    new_height = peak_height * envelope
                    banking_heights[j] = max(banking_heights[j], new_height)
                
                # End of if curve_length > 0
        
        # Handle curve that extends to end
        if in_curve:
            curve_end_idx = len(curve_flags) - 1
            
            # Clamp curve boundaries to transition points FIRST to prevent banking from crossing boundaries
            clamped_start_idx = curve_start_idx
            clamped_end_idx = curve_end_idx
            
            # Find last transition point before or at start
            for tp in sorted(transition_points, reverse=True):
                if tp < curve_start_idx:
                    clamped_start_idx = tp + 1  # Start after transition point
                    break
                elif tp == curve_start_idx:
                    # Transition point exactly at start - start after it
                    clamped_start_idx = tp + 1
                    break
            
            # Find first transition point at or after end, or within the curve
            for tp in sorted(transition_points):
                if tp > curve_end_idx:
                    clamped_end_idx = tp - 1  # End before transition point
                    break
                elif tp >= curve_start_idx and tp <= curve_end_idx:
                    # Transition point inside curve - clamp end to before it
                    clamped_end_idx = tp - 1
                    break
            
            # Ensure valid range
            clamped_start_idx = max(0, clamped_start_idx)
            clamped_end_idx = min(len(curve_flags) - 1, clamped_end_idx)
            
            # Only apply banking if clamped range is valid
            if clamped_start_idx <= clamped_end_idx:
                # Use clamped boundaries for all calculations
                curve_start_dist = distances[clamped_start_idx]
                curve_end_dist = distances[clamped_end_idx]
                curve_length = curve_end_dist - curve_start_dist
                
                if curve_length > 0:
                    # Start with full height, then reduce based on length only
                    peak_height = outer_edge_raise
                    
                    # Reduce height for short curves: reduce by ratio per meter below threshold
                    # Formula: new_h = h * (1.0 - (threshold - length) * ratio) (clamped to >= 0)
                    if curve_length < MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION:
                        meters_below_threshold = MIN_CURVE_LENGTH_THRESHOLD_FOR_REDUCTION - curve_length
                        reduction_factor = max(0.0, 1.0 - meters_below_threshold * CURVE_LENGTH_REDUCTION_RATIO)
                        peak_height *= reduction_factor
                    
                    # Check for adjacent opposite-direction curves
                    current_seg_idx = idx_to_seg.get(curve_start_idx, -1)
                    seg_dir = curve_segments_info[current_seg_idx]['direction'] if current_seg_idx >= 0 else 0.0
                    
                    # Determine fade zones at boundaries
                    fade_zone = 5.0  # meters
                    transition_distance = 5.0  # meters - minimum gap for fade
                    fade_start = 0.0
                    
                    # Check previous segment
                    if current_seg_idx > 0:
                        prev_seg = curve_segments_info[current_seg_idx - 1]
                        gap = curve_start_dist - prev_seg['end_dist']
                        # Only fade if there's a straight section (gap >= transition_distance)
                        if gap >= transition_distance:
                            prev_dir = prev_seg['direction']
                            if (prev_dir > 0 and seg_dir < 0) or (prev_dir < 0 and seg_dir > 0):
                                fade_start = min(fade_zone, curve_length / 2.0)
                    
                    # Apply parabolic arc envelope to clamped curve segment
                    # Envelope is recalculated over the clamped range to ensure smooth transitions
                    was_clamped_start = clamped_start_idx != curve_start_idx
                    
                    for j in range(clamped_start_idx, clamped_end_idx + 1):
                        # Position within clamped curve: 0.0 (clamped start) to 1.0 (clamped end)
                        # Clamp t to [0, 1] to handle floating point precision issues
                        t = (distances[j] - curve_start_dist) / curve_length
                        t = max(0.0, min(1.0, t))  # Clamp to [0, 1]
                        
                        # Parabolic arc: y = 4 * x * (1 - x)
                        # This gives a smooth curve that is 0 at t=0 and t=1, peaking at t=0.5
                        envelope = 4.0 * t * (1.0 - t)
                        
                        # Apply fade at start if opposite direction before (only if boundary wasn't clamped)
                        # Clamped boundaries already fade to 0 naturally via the envelope
                        if not was_clamped_start and fade_start > 0 and t < fade_start / curve_length:
                            fade_t = t * curve_length / fade_start
                            fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)  # smoothstep
                            envelope *= fade_factor
                        
                        # Blend with existing banking if curves overlap
                        new_height = peak_height * envelope
                        banking_heights[j] = max(banking_heights[j], new_height)
        
        # Smooth transitions between closely spaced curves
        # Find curve segments with their directions
        curve_segments = []
        in_curve = False
        seg_start_idx = None
        
        for i in range(len(curve_flags)):
            if curve_flags[i] and not in_curve:
                in_curve = True
                seg_start_idx = i
            elif not curve_flags[i] and in_curve:
                in_curve = False
                seg_end_idx = i - 1
                # Get representative direction for this segment (use middle or most common)
                mid_idx = (seg_start_idx + seg_end_idx) // 2
                seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
                curve_segments.append({
                    'start_idx': seg_start_idx,
                    'end_idx': seg_end_idx,
                    'start_dist': distances[seg_start_idx],
                    'end_dist': distances[seg_end_idx],
                    'direction': seg_dir
                })
        
        if in_curve:
            seg_end_idx = len(curve_flags) - 1
            mid_idx = (seg_start_idx + seg_end_idx) // 2
            seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
            curve_segments.append({
                'start_idx': seg_start_idx,
                'end_idx': seg_end_idx,
                'start_dist': distances[seg_start_idx],
                'end_dist': distances[seg_end_idx],
                'direction': seg_dir
            })
        
        # Handle transitions between curves - fade to zero when directions are opposite
        transition_distance = 5.0
        if len(curve_segments) > 1:
            for i in range(len(curve_segments) - 1):
                seg1 = curve_segments[i]
                seg2 = curve_segments[i + 1]
                gap = seg2['start_dist'] - seg1['end_dist']
                
                # Check if directions are opposite (left to right or right to left)
                dir1 = seg1['direction']
                dir2 = seg2['direction']
                opposite_directions = (dir1 > 0 and dir2 < 0) or (dir1 < 0 and dir2 > 0)
                
                if opposite_directions:
                    # Only fade banking to zero if there's a straight section between curves
                    # For back-to-back curves (gap < transition_distance), keep full banking
                    if gap >= transition_distance:
                        # Fade banking to zero at transition between opposite curves
                        # Use transition zone that extends into both segments
                        fade_zone = transition_distance
                    
                    for j in range(len(distances)):
                        if j >= len(banking_heights):
                            continue
                        dist = distances[j]
                        
                        # Fade out at end of seg1
                        if seg1['end_dist'] - fade_zone <= dist <= seg1['end_dist']:
                            if seg1['start_idx'] <= j <= seg1['end_idx']:
                                dist_from_end = seg1['end_dist'] - dist
                                fade_t = dist_from_end / fade_zone
                                fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)  # smoothstep
                                banking_heights[j] *= fade_factor
                        
                        # Zero banking in gap
                        elif seg1['end_dist'] < dist < seg2['start_dist']:
                            banking_heights[j] = 0.0
                        
                        # Fade in at start of seg2
                        elif seg2['start_dist'] <= dist <= seg2['start_dist'] + fade_zone:
                            if seg2['start_idx'] <= j <= seg2['end_idx']:
                                dist_from_start = dist - seg2['start_dist']
                                fade_t = dist_from_start / fade_zone
                                fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)  # smoothstep
                                banking_heights[j] *= fade_factor
                    # else: gap < transition_distance means back-to-back curves, keep full banking
                elif gap > 0 and gap < transition_distance:
                    # Same direction: blend normally
                    for j in range(len(distances)):
                        dist = distances[j]
                        if seg1['end_dist'] <= dist <= seg2['start_dist']:
                            t = (dist - seg1['end_dist']) / gap
                            smooth_t = t * t * (3.0 - 2.0 * t)  # smoothstep
                            h1 = banking_heights[seg1['end_idx']] if seg1['end_idx'] < len(banking_heights) else 0.0
                            h2 = banking_heights[seg2['start_idx']] if seg2['start_idx'] < len(banking_heights) else 0.0
                            blended = h1 * (1.0 - smooth_t) + h2 * smooth_t
                            banking_heights[j] = max(banking_heights[j], blended)
    
        return banking_heights
    
    def apply_curve_banking(self, geometry, dem_provider, z_raise, 
                           curve_window_distance, curve_threshold_angle,
                           outer_edge_raise, use_shapes=False, shape_distances_raw=None,
                           shape_heights=None, shape_table_lengths=None, sigmoid_lengths=None,
                           sigmoid_steepness_k=4.5, first_is_zero=False, last_is_end=False, feedback=None):
        """
        Apply curve detection and banking to a linestring.
        Detects curves by summing bisector angles over a distance window.
        For detected curves, raises the outer edge with a smooth parabolic arc.
        """
        # Validate geometry
        if geometry is None or geometry.isEmpty():
            return (geometry, None, None, None, None, None, None, None, [])
        
        # Store original geometry reference - we'll use it directly for XY positions
        original_geometry = geometry
        
        # Densify adaptively: 0.25m for curves, 1.25m for straights
        # This also performs curve detection on evenly-spaced 0.25m vertices
        # Returns original geometry reference + distances along it + curve flags
        curve_flags = []
        distances = []
        
        # Always use adaptive densification with fixed distance
        curve_threshold_angle_rad = math.radians(CURVE_THRESHOLD_ANGLE) if CURVE_THRESHOLD_ANGLE > 0 else 0.0
        original_geometry, curve_flags, curve_directions, distances, transition_points = self._densify_adaptive(geometry, curve_window_distance, curve_threshold_angle_rad, feedback)
        
        # Get vertices for DEM sampling - use distances to interpolate from original geometry
        # We need vertices for DEM sampling, but we'll use original geometry for XY in export
        num_vertices = len(distances) if distances else 0
        if num_vertices < 3:
            return (original_geometry, None, None, None, None, None, None, None, None)
        
        # Create temporary vertices list for DEM sampling (using original geometry)
        vertices = []
        for dist in distances:
            point = original_geometry.interpolate(dist)
            if point and not point.isEmpty():
                vertices.append(point.asPoint())
        
        if len(vertices) < 3:
            return (original_geometry, None, None, None, None, None, None, None, None)
        
        # Recalculate distances if not provided (shouldn't happen, but safety check)
        if not distances or len(distances) != len(vertices):
            distances = [0.0]
            total_length = 0.0
            for i in range(1, len(vertices)):
                prev_vertex = vertices[i-1]
                curr_vertex = vertices[i]
                dx = curr_vertex.x() - prev_vertex.x()
                dy = curr_vertex.y() - prev_vertex.y()
                segment_length = math.sqrt(dx*dx + dy*dy)
                total_length += segment_length
                distances.append(total_length)
        
        # Ensure curve_flags and curve_directions match vertex count
        if not curve_flags or len(curve_flags) != len(vertices):
            curve_flags = [False] * len(vertices)
        if not curve_directions or len(curve_directions) != len(vertices):
            curve_directions = [0.0] * len(vertices)
        
        # Step 1: Drape Z values from DEM at all vertices
        base_z_values = []
        for i, vertex in enumerate(vertices):
            x = vertex.x()
            y = vertex.y()
            val, ok = dem_provider.sample(QgsPointXY(x, y), 1)
            if ok and val is not None:
                base_z_values.append(float(val))
            else:
                base_z_values.append(0.0)
        
        # Step 2: Smooth the base Z values using 10 meter distance-based windows
        smoothing_distance = 10.0
        half_window = smoothing_distance / 2.0
        smoothed_z_values = []
        
        num_vertices = len(base_z_values)
        window_start_idx = 0
        
        for i in range(num_vertices):
            current_distance = distances[i]
            window_start_dist = current_distance - half_window
            window_end_dist = current_distance + half_window
            
            while window_start_idx < num_vertices and distances[window_start_idx] < window_start_dist:
                window_start_idx += 1
            
            window_z_values = []
            j = window_start_idx
            while j < num_vertices and distances[j] <= window_end_dist:
                window_z_values.append(base_z_values[j])
                j += 1
            
            if window_z_values:
                smoothed_z = sum(window_z_values) / len(window_z_values)
            else:
                smoothed_z = base_z_values[i]
            
            smoothed_z_values.append(smoothed_z)
        
        # Step 3: Curve directions are already calculated in _densify_adaptive
        # No need to recalculate - use the ones returned from adaptive densification
        # angle_sums are no longer needed (peak validation was removed)
        
        # Step 4: Apply shapes if provided
        shape_z_offsets = [0.0] * len(vertices)
        if use_shapes and shape_distances_raw and shape_heights and shape_table_lengths and sigmoid_lengths:
            total_length = distances[-1] if distances else 0.0
            shapes = list(zip([float(d) for d in shape_distances_raw], shape_heights, shape_table_lengths, sigmoid_lengths, range(len(shape_distances_raw))))
            
            for i, vertex in enumerate(vertices):
                distance_along = distances[i]
                shape_z = 0.0
                
                # Special case: First distance is 0
                if first_is_zero and len(shapes) > 0:
                    first_height = shape_heights[0]
                    first_table_len = shape_table_lengths[0]
                    first_sigmoid_len = sigmoid_lengths[0]
                    half_table_length = first_table_len / 2.0
                    half_wave_length = half_table_length + first_sigmoid_len
                    if distance_along <= half_wave_length:
                        if distance_along <= half_table_length:
                            shape_z_offset = first_height
                        else:
                            t = (distance_along - half_table_length) / first_sigmoid_len
                            if t <= 1.0:
                                sigmoid_val = self.normalized_sigmoid(t, sigmoid_steepness_k)
                                shape_z_offset = first_height * (1.0 - sigmoid_val)
                            else:
                                shape_z_offset = 0.0
                        shape_z_offset = max(0.0, shape_z_offset)
                        shape_z += shape_z_offset
                
                # Special case: Last distance is "end"
                if last_is_end and len(shapes) > 0:
                    last_height = shape_heights[-1]
                    last_table_len = shape_table_lengths[-1]
                    last_sigmoid_len = sigmoid_lengths[-1]
                    dist_from_end = total_length - distance_along
                    half_table_length = last_table_len / 2.0
                    half_wave_length = half_table_length + last_sigmoid_len
                    if dist_from_end >= 0 and dist_from_end <= half_wave_length:
                        if dist_from_end <= half_table_length:
                            shape_z_offset = last_height
                        elif dist_from_end <= half_table_length + last_sigmoid_len:
                            t = (dist_from_end - half_table_length) / last_sigmoid_len
                            if t <= 1.0:
                                sigmoid_val = self.normalized_sigmoid(t, sigmoid_steepness_k)
                                shape_z_offset = last_height * (1.0 - sigmoid_val)
                            else:
                                shape_z_offset = 0.0
                        else:
                            shape_z_offset = 0.0
                        shape_z_offset = max(0.0, shape_z_offset)
                        shape_z += shape_z_offset
                
                # Normal shapes
                for shape_dist, shape_height, table_len, sigmoid_len, orig_idx in shapes:
                    if first_is_zero and abs(shape_dist) < 1e-9:
                        continue
                    if last_is_end and abs(shape_dist - total_length) < 1e-9:
                        continue
                    
                    total_wave_length = sigmoid_len + table_len + sigmoid_len
                    dist_from_center = distance_along - shape_dist
                    
                    if abs(dist_from_center) <= total_wave_length / 2.0:
                        pos_from_wave_start = dist_from_center + total_wave_length / 2.0
                        
                        if pos_from_wave_start <= sigmoid_len:
                            t = pos_from_wave_start / sigmoid_len
                            if t <= 1.0:
                                sigmoid_val = self.normalized_sigmoid(t, sigmoid_steepness_k)
                                shape_z_offset = shape_height * sigmoid_val
                            else:
                                shape_z_offset = shape_height
                        elif pos_from_wave_start <= sigmoid_len + table_len:
                            shape_z_offset = shape_height
                        else:
                            t = (pos_from_wave_start - sigmoid_len - table_len) / sigmoid_len
                            if t <= 1.0:
                                sigmoid_val = self.normalized_sigmoid(t, sigmoid_steepness_k)
                                shape_z_offset = shape_height * (1.0 - sigmoid_val)
                            else:
                                shape_z_offset = 0.0
                        
                        shape_z_offset = max(0.0, shape_z_offset)
                        shape_z += shape_z_offset
                
                shape_z_offsets[i] = shape_z
        
        # Step 5: Create new points with raise, shapes, and curve banking applied
        # Store Z values in intermediate geometry for efficiency (includes shapes)
        # Original input geometry Z was stripped, but intermediate geometry Z is OK to use
        new_points = []
        for i, vertex in enumerate(vertices):
            x = vertex.x()  # XY from geometry (original input Z was stripped)
            y = vertex.y()  # XY from geometry (original input Z was stripped)
            base_z = smoothed_z_values[i]  # Base Z from DEM sampling
            
            # Base Z with raise and shapes - store in intermediate geometry for efficiency
            z = base_z + z_raise + shape_z_offsets[i]
            
            # Apply outer edge raise for curves (will be handled in ribbon generation)
            # For now, just store the base raised Z with shapes
            new_points.append(QgsPoint(x, y, z))
        
        # Recreate geometry (Z values stored but will be ignored in export)
        new_geometry = QgsGeometry.fromPolyline(new_points)
        
        # Create smooth parabolic arc envelope for banking
        banking_heights = self.create_parabolic_arc_envelope(
            curve_flags, distances, transition_points, 
            outer_edge_raise, curve_directions, feedback
        )
        
        if feedback:
            max_banking = max(banking_heights) if banking_heights else 0.0
            num_with_banking = sum(1 for h in banking_heights if h > 0.0)
            feedback.pushInfo(f"Banking heights calculated: max={max_banking:.3f}m, {num_with_banking} vertices with banking out of {len(banking_heights)} total")
        
        return (original_geometry, base_z_values, smoothed_z_values, curve_flags, 
                curve_directions, banking_heights, distances, shape_z_offsets)
    
    def _densify_adaptive(self, geometry, curve_window_distance, curve_threshold_angle, feedback=None):
        """
        Densify geometry adaptively: 0.25m for curves, 1.25m (5x) for straights.
        Process: 1) Densify to fixed 0.25m, 2) Calculate signed turning angles (+180/-180),
        3) Find transition points (direction changes), 4) Detect curves with single window (40% forward, 60% backward)
        clamped to transition points, 5) De-densify straights (keep every 5th vertex).
        Z values are NOT preserved - they will be sampled from DEM afterwards.
        
        Args:
            geometry: QgsGeometry linestring (2D or 3D, but only X/Y spacing is ensured)
            curve_window_distance: Distance window for curve detection (meters)
            curve_threshold_angle: Single threshold angle for curve detection (radians)
            feedback: Optional QgsProcessingFeedback for progress
            
        Returns:
            Tuple: (QgsGeometry, curve_flags list, curve_directions list, distances list, transition_points set)
        """
        # Fixed interval for curve detection (ensures even spacing)
        CURVE_DETECTION_INTERVAL = DENSIFY_DISTANCE  # meters
        STRAIGHT_DEDENSIFY_FACTOR = 5  # Keep every 5th vertex in straights
        
        if geometry is None or geometry.isEmpty():
            return (geometry, [], [], [], set())
        
        try:
            # Pass 1: Densify everything to fixed 0.25m for curve detection
            geometry_025m = self._densify_at_exact_intervals(geometry, CURVE_DETECTION_INTERVAL, feedback)
            
            # Get vertices at fixed spacing
            vertices_025m = list(geometry_025m.vertices())
            if len(vertices_025m) < 3:
                # Too few vertices, return as-is
                return (geometry, [False] * len(vertices_025m), [0.0] * len(vertices_025m), [0.0] * len(vertices_025m), set())
            
            # Use exact 0.25m intervals (not recalculated Euclidean distances)
            # This ensures consistent spacing for curve detection
            # Vertices are created at 0, DENSIFY_DISTANCE, 2*DENSIFY_DISTANCE, ..., total_length
            total_length = geometry.length()
            distances_025m = []
            for i in range(len(vertices_025m)):
                dist = min(i * CURVE_DETECTION_INTERVAL, total_length)
                distances_025m.append(dist)
            
            # Pass 2: Calculate signed turning angles (+180/-180 normalized) for each vertex
            # This is a vertex attribute: positive = left turn, negative = right turn, 0 = straight
            turning_angles = [0.0] * len(vertices_025m)  # Signed angles in radians
            
            for i in range(1, len(vertices_025m) - 1):
                    v_prev = vertices_025m[i-1]
                    v_curr = vertices_025m[i]
                    v_next = vertices_025m[i+1]
                    
                # Calculate vectors
                    dx1 = v_curr.x() - v_prev.x()
                    dy1 = v_curr.y() - v_prev.y()
                    dx2 = v_next.x() - v_curr.x()
                    dy2 = v_next.y() - v_curr.y()
                    
                    len1 = math.sqrt(dx1*dx1 + dy1*dy1)
                    len2 = math.sqrt(dx2*dx2 + dy2*dy2)
                    
                    if len1 < 1e-9 or len2 < 1e-9:
                        turning_angles[i] = 0.0
                        continue
                    
                # Normalize vectors
                    dx1 /= len1
                    dy1 /= len1
                    dx2 /= len2
                    dy2 /= len2
                    
                # Calculate angle between segments
                    dot = dx1*dx2 + dy1*dy2
                    dot = max(-1.0, min(1.0, dot))
                    angle = math.acos(dot)
                    
                    # Determine sign: left turn = positive, right turn = negative
                    cross = dx1*dy2 - dy1*dx2
                    sign = 1.0 if cross > 0 else (-1.0 if cross < 0 else 0.0)
                    
                    # Store signed turning angle (radians)
                    turning_angles[i] = angle * sign
            
            if feedback:
                min_angle = min(turning_angles) if turning_angles else 0.0
                max_angle = max(turning_angles) if turning_angles else 0.0
                feedback.pushInfo(f"Pass 2: Calculated signed turning angles for {len(turning_angles)} vertices")
                feedback.pushInfo(f"  Min angle: {min_angle:.4f} rad ({math.degrees(min_angle):.2f}°), Max angle: {max_angle:.4f} rad ({math.degrees(max_angle):.2f}°)")
            
            # Pass 3: Find transition points (where direction changes sign)
            # Transition points mark boundaries between opposite-direction curves
            transition_points = set()
            prev_sign = 0
            zero_run_start = None
            zero_run_last = None
            
            for i in range(1, len(turning_angles) - 1):
                sign = 1 if turning_angles[i] > 0 else (-1 if turning_angles[i] < 0 else 0)
                
                if sign == 0:
                    # Track run of zeros
                    if zero_run_start is None:
                        zero_run_start = i
                    zero_run_last = i  # Update last zero in run
                else:
                    # Non-zero sign
                    if zero_run_start is not None:
                        # Multiple zeros in a row - use LAST zero as transition point
                        if zero_run_last is not None:
                            transition_points.add(zero_run_last)
                        zero_run_start = None
                        zero_run_last = None
                    
                    # Check for sign change
                    if prev_sign != 0 and sign != prev_sign:
                        # Direction changed - this is a transition point
                        transition_points.add(i)
                    
                    prev_sign = sign
            
            # Handle final zero run if exists
            if zero_run_start is not None and zero_run_last is not None:
                transition_points.add(zero_run_last)
            
            # Convert to sorted list for easier lookup
            transition_points_list = sorted(transition_points)
            
            if feedback:
                feedback.pushInfo(f"Pass 3: Found {len(transition_points_list)} transition points")
                if transition_points_list:
                    tp_info = []
                    for tp in transition_points_list[:20]:
                        tp_dist = distances_025m[tp] if tp < len(distances_025m) else 0.0
                        tp_info.append(f"v{tp}({tp_dist:.2f}m)")
                    feedback.pushInfo(f"  Transition points: {', '.join(tp_info)}{'...' if len(transition_points_list) > 20 else ''}")
            
            # Pass 4: Curve detection with single window (40% forward, 60% backward)
            # Clamp curve boundaries to transition points
            curve_flags_025m = [False] * len(vertices_025m)
            curve_directions_025m = [0.0] * len(vertices_025m)
            
            if curve_window_distance > 0 and curve_threshold_angle > 0:
                forward_dist = CURVE_WINDOW_FORWARD_RATIO * curve_window_distance
                backward_dist = CURVE_WINDOW_BACKWARD_RATIO * curve_window_distance
                
                if feedback:
                    feedback.pushInfo(f"Pass 4: Curve detection with single window")
                    feedback.pushInfo(f"  Window: {forward_dist:.2f}m forward ({CURVE_WINDOW_FORWARD_RATIO*100:.0f}%), {backward_dist:.2f}m backward ({CURVE_WINDOW_BACKWARD_RATIO*100:.0f}%), total {curve_window_distance}m")
                    feedback.pushInfo(f"  Threshold: {math.degrees(curve_threshold_angle):.2f}° ({curve_threshold_angle:.4f} rad)")
                
                # Calculate window sums for each vertex using sliding window (more efficient)
                window_sums = [0.0] * len(vertices_025m)
                
                for i in range(len(vertices_025m)):
                    current_distance = distances_025m[i]
                    
                    # Window: 40% forward, 60% backward
                    window_start_dist = max(0.0, current_distance - backward_dist)
                    window_end_dist = min(total_length, current_distance + forward_dist)
                    
                    # Find start and end indices for window (more efficient than checking all vertices)
                    window_start_idx = 0
                    window_end_idx = len(vertices_025m) - 1
                    
                    # Binary search would be better, but linear search is simpler and fast enough for 0.25m spacing
                    for j in range(len(vertices_025m)):
                        if distances_025m[j] >= window_start_dist:
                            window_start_idx = j
                            break
                    
                    for j in range(len(vertices_025m) - 1, -1, -1):
                        if distances_025m[j] <= window_end_dist:
                            window_end_idx = j
                            break
                    
                    # Sum signed turning angles in window
                    angle_sum = 0.0
                    for j in range(window_start_idx, window_end_idx + 1):
                        if 0 < j < len(vertices_025m) - 1:
                            angle_sum += turning_angles[j]
                    
                    window_sums[i] = angle_sum
                
                # Detect curves: find where threshold is exceeded, then clamp to transition points
                # Group consecutive vertices that exceed threshold into curve segments
                in_potential_curve = False
                potential_curve_start = None
                curve_count = 0
                
                for i in range(len(vertices_025m)):
                    if abs(window_sums[i]) >= curve_threshold_angle:
                        if not in_potential_curve:
                            # Start of potential curve
                            in_potential_curve = True
                            potential_curve_start = i
                    else:
                        if in_potential_curve:
                            # End of potential curve - process it
                            potential_curve_end = i - 1
                            
                            # Find peak (maximum |sum|) within potential curve
                            peak_idx = potential_curve_start
                            peak_sum = abs(window_sums[potential_curve_start])
                            for j in range(potential_curve_start, potential_curve_end + 1):
                                if abs(window_sums[j]) > peak_sum:
                                    peak_sum = abs(window_sums[j])
                                    peak_idx = j
                            
                            # Determine direction from peak
                            direction = 1.0 if window_sums[peak_idx] > 0 else -1.0
                            
                            # Find nearest transition points to clamp boundaries
                            # Clamp to prevent curve from crossing any transition points
                            # BUT ensure peak is always within clamped range
                            clamped_start = potential_curve_start
                            clamped_end = potential_curve_end
                            
                            # Find last transition point before peak (but not before potential_curve_start)
                            # This ensures we don't clamp too far back
                            for tp in reversed(transition_points_list):
                                if tp < peak_idx and tp >= potential_curve_start:
                                    # Transition point between start and peak - start after it
                                    clamped_start = tp + 1
                                    break
                                elif tp < potential_curve_start:
                                    # Transition point before potential start - don't clamp to it
                                    # (keep original start)
                                    break
                            
                            # Find first transition point after peak (but not after potential_curve_end)
                            # This ensures we don't clamp too far forward
                            for tp in transition_points_list:
                                if tp > peak_idx and tp <= potential_curve_end:
                                    # Transition point between peak and end - end before it
                                    clamped_end = tp - 1
                                    break
                                elif tp > potential_curve_end:
                                    # Transition point after potential end - don't clamp to it
                                    # (keep original end)
                                    break
                            
                            # Ensure peak is always within clamped range
                            clamped_start = min(clamped_start, peak_idx)
                            clamped_end = max(clamped_end, peak_idx)
                            
                            # Ensure valid range
                            clamped_start = max(0, clamped_start)
                            clamped_end = min(len(vertices_025m) - 1, clamped_end)
                            
                            # Flag vertices in clamped range as curve
                            if clamped_start <= clamped_end:
                                for j in range(clamped_start, clamped_end + 1):
                                    curve_flags_025m[j] = True
                                    curve_directions_025m[j] = direction
                                
                                if feedback:
                                    dir_str = "LEFT" if direction > 0 else "RIGHT"
                                    start_dist = distances_025m[clamped_start] if clamped_start < len(distances_025m) else 0.0
                                    end_dist = distances_025m[clamped_end] if clamped_end < len(distances_025m) else 0.0
                                    peak_dist = distances_025m[peak_idx] if peak_idx < len(distances_025m) else 0.0
                                    pot_start_dist = distances_025m[potential_curve_start] if potential_curve_start < len(distances_025m) else 0.0
                                    pot_end_dist = distances_025m[potential_curve_end] if potential_curve_end < len(distances_025m) else 0.0
                                    was_clamped = (potential_curve_start != clamped_start) or (potential_curve_end != clamped_end)
                                    clamp_info = f" [CLAMPED from v{potential_curve_start}-{potential_curve_end} ({pot_start_dist:.2f}-{pot_end_dist:.2f}m)]" if was_clamped else ""
                                    feedback.pushInfo(f"  Curve {dir_str}: vertices {clamped_start}-{clamped_end} (peak at {peak_idx}), distances {start_dist:.2f}-{end_dist:.2f}m (peak {peak_dist:.2f}m), length {end_dist-start_dist:.2f}m{clamp_info}")
                                    curve_count += 1
                            
                            in_potential_curve = False
                
                # Handle curve that extends to end
                if in_potential_curve:
                    potential_curve_end = len(vertices_025m) - 1
                    
                    # Find peak
                    peak_idx = potential_curve_start
                    peak_sum = abs(window_sums[potential_curve_start])
                    for j in range(potential_curve_start, potential_curve_end + 1):
                        if abs(window_sums[j]) > peak_sum:
                            peak_sum = abs(window_sums[j])
                            peak_idx = j
                    
                    direction = 1.0 if window_sums[peak_idx] > 0 else -1.0
                    
                    # Clamp start to transition point (but ensure peak is included)
                    clamped_start = potential_curve_start
                    for tp in reversed(transition_points_list):
                        if tp < peak_idx and tp >= potential_curve_start:
                            # Transition point between start and peak - start after it
                            clamped_start = tp + 1
                            break
                        elif tp < potential_curve_start:
                            # Transition point before potential start - don't clamp to it
                            break
                    
                    # Clamp end to transition point if any exists within or after curve
                    clamped_end = potential_curve_end
                    for tp in transition_points_list:
                        if tp > peak_idx and tp <= potential_curve_end:
                            # Transition point between peak and end - end before it
                            clamped_end = tp - 1
                            break
                        elif tp > potential_curve_end:
                            # Transition point after potential end - don't clamp to it
                            break
                    
                    # Ensure peak is always within clamped range
                    clamped_start = min(clamped_start, peak_idx)
                    clamped_end = max(clamped_end, peak_idx)
                    
                    clamped_start = max(0, clamped_start)
                    
                    # Flag vertices
                    for j in range(clamped_start, potential_curve_end + 1):
                        curve_flags_025m[j] = True
                        curve_directions_025m[j] = direction
                    
                    if feedback:
                        dir_str = "LEFT" if direction > 0 else "RIGHT"
                        start_dist = distances_025m[clamped_start] if clamped_start < len(distances_025m) else 0.0
                        end_dist = distances_025m[potential_curve_end] if potential_curve_end < len(distances_025m) else 0.0
                        peak_dist = distances_025m[peak_idx] if peak_idx < len(distances_025m) else 0.0
                        pot_start_dist = distances_025m[potential_curve_start] if potential_curve_start < len(distances_025m) else 0.0
                        was_clamped = (potential_curve_start != clamped_start) or (potential_curve_end != clamped_end)
                        clamp_info = f" [CLAMPED from v{potential_curve_start}-{potential_curve_end} ({pot_start_dist:.2f}-{end_dist:.2f}m)]" if was_clamped else ""
                        feedback.pushInfo(f"  Curve {dir_str} (to end): vertices {clamped_start}-{potential_curve_end} (peak at {peak_idx}), distances {start_dist:.2f}-{end_dist:.2f}m (peak {peak_dist:.2f}m), length {end_dist-start_dist:.2f}m{clamp_info}")
                        curve_count += 1
                
                if feedback:
                    total_curve_vertices = sum(1 for f in curve_flags_025m if f)
                    feedback.pushInfo(f"Pass 4 Summary: Detected {curve_count} curve segments, {total_curve_vertices} curve vertices out of {len(vertices_025m)} total vertices")
            
            # Pass 5: De-densify straights (keep every 5th vertex, keep all curve vertices)
            # Instead of creating new geometry, we keep distances along ORIGINAL geometry
            keep_indices = []
            
            # First pass: identify curve segment boundaries (where curve flags start/end)
            # Only mark straight vertices adjacent to curve segments (curve vertices are already kept)
            curve_boundaries = set()
            in_curve = False
            for i in range(len(curve_flags_025m)):
                if curve_flags_025m[i] and not in_curve:
                    # Start of curve segment
                    in_curve = True
                    # Keep straight vertex before curve (if it exists and is not a curve)
                    if i > 0 and not curve_flags_025m[i - 1]:
                        curve_boundaries.add(i - 1)
                elif not curve_flags_025m[i] and in_curve:
                    # End of curve segment
                    in_curve = False
                    # Keep straight vertex after curve (already confirmed not a curve by elif condition)
                    curve_boundaries.add(i)
            
            # Second pass: decide which vertices to keep
            for i in range(len(vertices_025m)):
                keep_vertex = False
                
                # Always keep first and last vertices
                if i == 0 or i == len(vertices_025m) - 1:
                    keep_vertex = True
                # Keep all curve vertices
                elif curve_flags_025m[i]:
                    keep_vertex = True
                # Keep vertices at curve segment boundaries (for smooth transitions)
                elif i in curve_boundaries:
                    keep_vertex = True
                # For straights: keep every 5th vertex (0, 5, 10, 15, ...)
                elif i % STRAIGHT_DEDENSIFY_FACTOR == 0:
                    keep_vertex = True
                
                if keep_vertex:
                    keep_indices.append(i)
            
            # Map curve flags, directions, and distances to kept indices
            # IMPORTANT: distances_025m are distances along the 0.25m-densified geometry
            # We need to convert these to distances along the ORIGINAL geometry
            # Since geometry_025m was created by interpolating original geometry at exact distances,
            # distances_025m[i] corresponds to distance i * 0.25 along the original geometry
            curve_flags_adaptive = [curve_flags_025m[i] for i in keep_indices]
            curve_directions_adaptive = [curve_directions_025m[i] for i in keep_indices]
            distances_adaptive = [distances_025m[i] for i in keep_indices]
            
            # Map transition points to adaptive indices (find closest kept index for each transition point)
            transition_points_adaptive = set()
            for tp in transition_points_list:
                # Find closest kept index to this transition point
                closest_idx = 0
                min_dist = abs(keep_indices[0] - tp)
                for idx, kept_i in enumerate(keep_indices):
                    dist = abs(kept_i - tp)
                    if dist < min_dist:
                        min_dist = dist
                        closest_idx = idx
                transition_points_adaptive.add(closest_idx)
            
            # Return ORIGINAL geometry reference (not new geometry)
            # XY positions will be calculated directly from original geometry using distances
            if feedback:
                num_curves = sum(1 for f in curve_flags_adaptive if f)
                feedback.pushInfo(f"Adaptive densification: {len(vertices_025m)} vertices at 0.25m → {len(keep_indices)} distances (curves: {num_curves} at 0.25m, straights: {len(keep_indices) - num_curves} at 1.25m)")
                feedback.pushInfo(f"Found {len(transition_points_adaptive)} transition points")
                
                # Summary of detected curves
                in_curve = False
                curve_num = 0
                for i in range(len(curve_flags_adaptive)):
                    if curve_flags_adaptive[i] and not in_curve:
                        in_curve = True
                        curve_start = i
                        curve_dir = curve_directions_adaptive[i]
                    elif not curve_flags_adaptive[i] and in_curve:
                        in_curve = False
                        curve_end = i - 1
                        curve_num += 1
                        dir_str = "LEFT" if curve_dir > 0 else "RIGHT"
                        start_dist = distances_adaptive[curve_start] if curve_start < len(distances_adaptive) else 0.0
                        end_dist = distances_adaptive[curve_end] if curve_end < len(distances_adaptive) else 0.0
                        feedback.pushInfo(f"  Curve {curve_num} ({dir_str}): adaptive indices {curve_start}-{curve_end}, distances {start_dist:.2f}-{end_dist:.2f}m")
                
                if in_curve:
                    curve_num += 1
                    curve_end = len(curve_flags_adaptive) - 1
                    dir_str = "LEFT" if curve_dir > 0 else "RIGHT"
                    start_dist = distances_adaptive[curve_start] if curve_start < len(distances_adaptive) else 0.0
                    end_dist = distances_adaptive[curve_end] if curve_end < len(distances_adaptive) else 0.0
                    feedback.pushInfo(f"  Curve {curve_num} ({dir_str}): adaptive indices {curve_start}-{curve_end}, distances {start_dist:.2f}-{end_dist:.2f}m")
            
            return (geometry, curve_flags_adaptive, curve_directions_adaptive, distances_adaptive, transition_points_adaptive)
                
        except Exception as e:
            if feedback:
                feedback.pushWarning(f"Error in adaptive densification: {str(e)}, falling back to uniform 0.25m densification")
            # Fallback: return original geometry with distances at 0.25m intervals
            total_length = geometry.length()
            if total_length > 0:
                fallback_distances = []
                num_intervals = int(math.ceil(total_length / CURVE_DETECTION_INTERVAL))
                for i in range(num_intervals + 1):
                    dist = i * CURVE_DETECTION_INTERVAL
                    if dist <= total_length:
                        fallback_distances.append(dist)
                if fallback_distances[-1] < total_length - 1e-6:
                    fallback_distances.append(total_length)
                return (geometry, [False] * len(fallback_distances), [0.0] * len(fallback_distances), fallback_distances, set())
            return (geometry, [], [], [0.0], set())
    
    def _densify_at_exact_intervals(self, geometry, interval_distance, feedback=None):
        """
        Densify geometry to have vertices at exact distance intervals in X/Y only (2D).
        Ensures vertices are placed at 0, interval_distance, 2*interval_distance, etc.
        along the line, guaranteeing at least one vertex at each interval.
        Z values are NOT preserved - they will be sampled from DEM afterwards.
        
        Args:
            geometry: QgsGeometry linestring (2D or 3D, but only X/Y spacing is ensured)
            interval_distance: Desired distance between vertices (e.g., 0.25m)
            feedback: Optional QgsProcessingFeedback for progress
            
        Returns:
            QgsGeometry with vertices at exact intervals (2D, no Z values)
        """
        if geometry is None or geometry.isEmpty() or interval_distance <= 0:
            return geometry
        
        try:
            # Get total length (2D length in X/Y plane)
            total_length = geometry.length()
            if total_length <= 0:
                return geometry
            
            # Calculate number of intervals
            num_intervals = int(math.ceil(total_length / interval_distance))
            
            # Create list of exact distances: 0, interval, 2*interval, ..., total_length
            exact_distances = []
            for i in range(num_intervals + 1):
                dist = i * interval_distance
                if dist <= total_length:
                    exact_distances.append(dist)
            
            # Always include the end point
            if exact_distances[-1] < total_length - 1e-6:
                exact_distances.append(total_length)
            
            # Interpolate points at exact distances (2D only, no Z)
            new_points = []
            for dist in exact_distances:
                # Use QGIS interpolation to get point at exact distance (2D)
                point = geometry.interpolate(dist)
                if point and not point.isEmpty():
                    point_xy = point.asPoint()
                    # Create 2D point only (no Z) - Z will come from DEM sampling
                    new_points.append(QgsPoint(point_xy.x(), point_xy.y()))
            
            # Create new geometry from interpolated points (2D)
            if len(new_points) >= 2:
                new_geometry = QgsGeometry.fromPolyline(new_points)
                return new_geometry
            else:
                return geometry
                
        except Exception as e:
            if feedback:
                feedback.pushWarning(f"Error in exact interval densification: {str(e)}")
            # Fallback to standard densification
            try:
                return geometry.densifyByDistance(interval_distance)
            except Exception:
                return geometry
    
    def _resample_cross_section_positions(self, positions_list, target_count):
        """
        Resample a list of 3D positions to a target count using linear interpolation.
        This ensures smooth transitions when connecting cross-sections with different vertex counts.
        
        Args:
            positions_list: List of [x, y, z] positions
            target_count: Desired number of positions
            
        Returns:
            List of [x, y, z] positions with target_count elements
        """
        if len(positions_list) == target_count:
            return positions_list
        
        if len(positions_list) == 1:
            # Single position: repeat it
            return positions_list * target_count
        
        if target_count == 1:
            # Return middle position
            return [positions_list[len(positions_list) // 2]]
        
        # Linear interpolation along the cross-section
        # CRITICAL: Always preserve first and last vertices exactly to prevent gaps
        resampled = []
        for i in range(target_count):
            # Always use exact edge vertices
            if i == 0:
                resampled.append(positions_list[0])
                continue
            elif i == target_count - 1:
                resampled.append(positions_list[-1])
                continue
            
            # Interpolate intermediate vertices
            t = i / (target_count - 1) if target_count > 1 else 0.0
            source_idx = t * (len(positions_list) - 1)
            idx_low = int(source_idx)
            idx_high = min(idx_low + 1, len(positions_list) - 1)
            alpha = source_idx - idx_low
            
            # Interpolate 3D position
            p0 = positions_list[idx_low]
            p1 = positions_list[idx_high]
            resampled.append([
                p0[0] * (1 - alpha) + p1[0] * alpha,
                p0[1] * (1 - alpha) + p1[1] * alpha,
                p0[2] * (1 - alpha) + p1[2] * alpha
            ])
        return resampled
    
    def export_ribbon_mesh(self, feature, original_geometry, base_z_values, smoothed_z_values,
                          curve_flags, curve_directions, banking_heights, distances,
                          shape_z_offsets, ribbon_width, transversal_concavity, export_format,
                          project_home, texture_file=None, input_layer=None, dem_layer=None,
                          z_raise=1.0, feedback=None):
        """
        Export a ribbon mesh with curve banking and optional transversal concavity.
        Uses original geometry directly for XY positions to maintain precision.
        For curved sections, raises the outer edge using smooth parabolic arc envelope.
        Transversal concavity creates a concave (facing down) arc across the width.
        """
        if original_geometry is None or original_geometry.isEmpty():
            return
        
        if not distances or len(distances) == 0:
            return
        
        num_vertices = len(distances)
        if len(base_z_values) != num_vertices or len(smoothed_z_values) != num_vertices:
            if feedback:
                feedback.pushWarning(f"Feature {feature.id()}: Mismatch in vertex counts, skipping ribbon export.")
            return
        
        # Generate ribbon mesh vertices
        ribbon_vertices = []
        base_half_width = ribbon_width / 2.0
        bottom_offset = 1.0
        
        # Small delta for calculating tangent from original geometry
        TANGENT_DELTA = 0.01  # 1cm
        
        # Pre-calculate perpendicular vectors directly from original geometry
        # This ensures perfect alignment with the original line
        vertex_perpendiculars = []
        for i in range(num_vertices):
            distance = distances[i]
            
            # Get point on original geometry at this distance
            point_geom = original_geometry.interpolate(distance)
            if not point_geom or point_geom.isEmpty():
                # Fallback to default
                vertex_perpendiculars.append([1.0, 0.0])
                continue
            
            point = point_geom.asPoint()
            
            # Calculate tangent by getting point slightly ahead
            next_distance = min(distance + TANGENT_DELTA, distances[-1] if distances else distance + TANGENT_DELTA)
            if next_distance <= distance:
                # At end, use point slightly before
                prev_distance = max(distance - TANGENT_DELTA, 0.0)
                prev_point_geom = original_geometry.interpolate(prev_distance)
                if prev_point_geom and not prev_point_geom.isEmpty():
                    prev_point = prev_point_geom.asPoint()
                    dx = point.x() - prev_point.x()
                    dy = point.y() - prev_point.y()
                else:
                    vertex_perpendiculars.append([1.0, 0.0])
                    continue
            else:
                next_point_geom = original_geometry.interpolate(next_distance)
                if next_point_geom and not next_point_geom.isEmpty():
                    next_point = next_point_geom.asPoint()
                    dx = next_point.x() - point.x()
                    dy = next_point.y() - point.y()
                else:
                    vertex_perpendiculars.append([1.0, 0.0])
                    continue
            
            # Normalize tangent
            length = math.sqrt(dx*dx + dy*dy)
            if length < 1e-9:
                # Degenerate, use previous or default
                if len(vertex_perpendiculars) > 0:
                    vertex_perpendiculars.append(vertex_perpendiculars[-1])
                else:
                    vertex_perpendiculars.append([1.0, 0.0])
                continue
            
            dx /= length
            dy /= length
            
            # Perpendicular vector (left = counterclockwise 90 degrees)
            nx = -dy
            ny = dx
            vertex_perpendiculars.append([nx, ny])
        
        # Pre-calculate curve segment positions for concavity scaling
        # Track which curve segment each vertex belongs to and its position within that segment
        curve_segment_positions = [0.0] * num_vertices  # 0.0 = start/end of curve, 1.0 = middle
        transition_distance = 2.0  # Distance in meters for transition zones at curve boundaries
        in_curve = False
        curve_start_idx = None
        curve_segments = []  # Store curve segment info for transition zone processing
        
        for i in range(num_vertices):
            if curve_flags[i] and not in_curve:
                # Start of new curve segment
                in_curve = True
                curve_start_idx = i
            elif not curve_flags[i] and in_curve:
                # End of curve segment
                in_curve = False
                curve_end_idx = i - 1
                curve_start_dist = distances[curve_start_idx]
                curve_end_dist = distances[curve_end_idx]
                curve_length = curve_end_dist - curve_start_dist
                
                if curve_length > 0:
                    # Calculate position within curve segment for each vertex
                    # Use full concavity factor (1.0) throughout the curve for consistent application
                    # The transition zones will handle the fade-in/out at boundaries
                    for j in range(curve_start_idx, curve_end_idx + 1):
                        # Use full concavity factor in curves (1.0)
                        # Transition zones in straight sections will fade from 0 to 1.0 at entrance
                        # and from 1.0 to 0 at exit
                        curve_segment_positions[j] = 1.0
                    
                    # Store curve segment info for transition zone processing
                    curve_segments.append({
                        'start_idx': curve_start_idx,
                        'end_idx': curve_end_idx,
                        'start_dist': curve_start_dist,
                        'end_dist': curve_end_dist
                    })
        
        # Handle curve that extends to end
        if in_curve:
            curve_end_idx = len(curve_flags) - 1
            curve_start_dist = distances[curve_start_idx]
            curve_end_dist = distances[curve_end_idx]
            curve_length = curve_end_dist - curve_start_dist
            
            if curve_length > 0:
                # Use full concavity factor (1.0) throughout the curve
                for j in range(curve_start_idx, curve_end_idx + 1):
                    curve_segment_positions[j] = 1.0
                
                # Store curve segment info
                curve_segments.append({
                    'start_idx': curve_start_idx,
                    'end_idx': curve_end_idx,
                    'start_dist': curve_start_dist,
                    'end_dist': curve_end_dist
                })
        
        # Apply transition zones: fade concavity only in straight sections between curves
        # If curves are back-to-back (no straight), keep concavity at 1.0
        for seg_idx, seg in enumerate(curve_segments):
            curve_start_dist = seg['start_dist']
            curve_end_dist = seg['end_dist']
            
            # Check if there's a previous curve segment
            has_prev_curve = seg_idx > 0
            if has_prev_curve:
                prev_seg = curve_segments[seg_idx - 1]
                gap_before = curve_start_dist - prev_seg['end_dist']
            else:
                gap_before = float('inf')  # No previous curve, treat as large gap
            
            # Check if there's a next curve segment
            has_next_curve = seg_idx < len(curve_segments) - 1
            if has_next_curve:
                next_seg = curve_segments[seg_idx + 1]
                gap_after = next_seg['start_dist'] - curve_end_dist
            else:
                gap_after = float('inf')  # No next curve, treat as large gap
            
            for j in range(num_vertices):
                dist = distances[j]
                
                # Skip if already in a curve (keep curve's own concavity profile)
                if curve_segment_positions[j] > 0.0:
                    continue
                
                # Before curve start: fade concavity from 0 to 1.0 only if there's a straight section
                if gap_before >= transition_distance:
                    if curve_start_dist - transition_distance <= dist < curve_start_dist:
                        # Distance from curve start (positive = before curve)
                        dist_from_curve = curve_start_dist - dist
                        # Fade factor: 0.0 at transition_distance away, 1.0 at curve boundary
                        fade_t = dist_from_curve / transition_distance
                        # Smoothstep for smooth fade-in
                        fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)
                        # Apply concavity factor that fades from 0 to 1.0
                        curve_segment_positions[j] = fade_factor
                
                # After curve end: fade concavity from 1.0 to 0 only if there's a straight section
                if gap_after >= transition_distance:
                    if curve_end_dist < dist <= curve_end_dist + transition_distance:
                        # Distance from curve end (positive = after curve)
                        dist_from_curve = dist - curve_end_dist
                        # Fade factor: 1.0 at curve boundary, 0.0 at transition_distance away
                        fade_t = 1.0 - (dist_from_curve / transition_distance)
                        # Smoothstep for smooth fade-out
                        fade_factor = fade_t * fade_t * (3.0 - 2.0 * fade_t)
                        # Apply concavity factor that fades from 1.0 to 0
                        curve_segment_positions[j] = fade_factor
        
        # Generate ribbon vertices with curve banking
        # Use original geometry directly for XY positions to maintain precision
        for i in range(num_vertices):
            distance = distances[i]
            
            # Get exact XY position from original geometry
            point_geom = original_geometry.interpolate(distance)
            if not point_geom or point_geom.isEmpty():
                continue
            
            point = point_geom.asPoint()
            x = point.x()  # Exact XY from original geometry
            y = point.y()  # Exact XY from original geometry
            base_z = smoothed_z_values[i]  # Base elevation from DEM sampling
            
            # Calculate upper Z from DEM-sampled base + raise + shapes
            shape_z = shape_z_offsets[i] if i < len(shape_z_offsets) else 0.0
            upper_z = base_z + z_raise + shape_z
            
            # Use the perpendicular assigned to this vertex
            nx, ny = vertex_perpendiculars[i]
            
            # Determine if this vertex is in a curve and which side is outer
            is_curve = curve_flags[i] if i < len(curve_flags) else False
            curve_dir = curve_directions[i] if i < len(curve_directions) else 0.0
            
            # Use base width (no widening)
            half_width = base_half_width
            
            # Calculate lower Z
            lower_z = base_z - bottom_offset
            
            # Use pre-calculated smooth parabolic arc banking height
            banking_height = banking_heights[i] if i < len(banking_heights) else 0.0
            
            # Determine which side is outer for curves and apply banking FIRST
            # curve_dir: 1.0 = left curve (turning left, outer edge is RIGHT), -1.0 = right curve (turning right, outer edge is LEFT)
            if banking_height > 0.0 and curve_dir > 0:
                # Left curve: turning left, so outer edge is on the RIGHT
                # Raise right side proportionally
                upper_left_z = upper_z
                upper_right_z = upper_z + banking_height
            elif banking_height > 0.0 and curve_dir < 0:
                # Right curve: turning right, so outer edge is on the LEFT
                # Raise left side proportionally
                upper_left_z = upper_z + banking_height
                upper_right_z = upper_z
            else:
                # Straight or below threshold: no banking
                upper_left_z = upper_z
                upper_right_z = upper_z
            
            # Calculate sidewall heights AFTER banking is applied, so back slopes are always 45°
            sidewall_height_left = upper_left_z - lower_z
            sidewall_height_right = upper_right_z - lower_z
            horizontal_offset_left = sidewall_height_left  # 45° incline
            horizontal_offset_right = sidewall_height_right  # 45° incline
            
            # Create vertices for this cross-section
            # Note: +nx points left, -nx points right (matching original script and non-concavity case)
            lower_left = [x + nx * (half_width + horizontal_offset_left), y + ny * (half_width + horizontal_offset_left), lower_z]
            lower_right = [x - nx * (half_width + horizontal_offset_right), y - ny * (half_width + horizontal_offset_right), lower_z]
            
            # Determine inner lower edge Z (the lower edge on the inner side of the curve)
            if banking_height > 0.0 and curve_dir > 0:
                # Left curve: inner edge is LEFT
                inner_lower_z = lower_left[2]
            elif banking_height > 0.0 and curve_dir < 0:
                # Right curve: inner edge is RIGHT
                inner_lower_z = lower_right[2]
            else:
                # Straight: use either (they're the same)
                inner_lower_z = lower_left[2]
            
            # Apply transversal concavity with smooth transitions at curve boundaries
            # Use curve_segment_positions which includes transition zones (fades to 0 at boundaries)
            if transversal_concavity > 0.0 and curve_segment_positions[i] > 0.0:
                # Calculate effective concavity: scaled by banking height and curve segment position
                # curve_segment_positions[i] is 0.0 at start/end, peaks at 1.0 in middle (quadratic envelope)
                # Scale by banking_height so concavity increases proportionally with banking
                # Use a reference banking height (e.g., 1.0m) to normalize the scaling
                reference_banking = 1.0  # Reference banking height for scaling
                banking_scale = min(1.0, banking_height / reference_banking) if reference_banking > 0 else 0.0
                effective_concavity = transversal_concavity * curve_segment_positions[i] * banking_scale
                
                # Create multiple upper vertices across width with concave profile
                num_transversal_segments = 5  # 5 segments = 6 vertices (left + 4 intermediate + right)
                upper_vertices = []
                
                for j in range(num_transversal_segments + 1):
                    # Position across width: 0.0 (left edge) to 1.0 (right edge)
                    t = j / num_transversal_segments
                    
                    # Interpolate edge Z values (accounting for banking)
                    # At t=0: left edge, at t=1: right edge
                    edge_z = upper_left_z * (1.0 - t) + upper_right_z * t
                    
                    # Concave parabolic profile: center is lower (facing down)
                    # Formula: z = edge_z - concavity * (4 * t * (1 - t))
                    # At edges (t=0, t=1): multiplier = 0, so z = edge_z (higher)
                    # At center (t=0.5): multiplier = 1.0, so z = edge_z - concavity (lower)
                    concavity_factor = 4.0 * t * (1.0 - t)
                    z = edge_z - effective_concavity * concavity_factor
                    
                    # Ensure concavity never goes below inner lower edge Z
                    z = max(z, inner_lower_z)
                    
                    # Calculate position across width
                    # From +half_width (left, t=0) to -half_width (right, t=1) to match non-concavity case
                    # Non-concavity: left uses +nx, right uses -nx
                    offset = half_width - (2.0 * half_width * t)
                    vx = x + nx * offset
                    vy = y + ny * offset
                    
                    upper_vertices.append([vx, vy, z])
                
                ribbon_vertices.append({
                    'lower_left': lower_left,
                    'lower_right': lower_right,
                    'upper_vertices': upper_vertices,
                    'num_upper_vertices': len(upper_vertices),
                    'has_concavity': True
                })
            else:
                # No concavity: use simple 2-vertex approach
                upper_left = [x + nx * half_width, y + ny * half_width, upper_left_z]
                upper_right = [x - nx * half_width, y - ny * half_width, upper_right_z]
                
                ribbon_vertices.append({
                    'lower_left': lower_left,
                    'lower_right': lower_right,
                    'upper_left': upper_left,
                    'upper_right': upper_right,
                    'has_concavity': False
                })
        
        if len(ribbon_vertices) < 2:
            return
        
        # Build mesh: triangulate ribbon with all faces
        positions = []
        indices = []
        texcoords = []
        
        # Helper function to convert QGIS coords to GLTF coords
        def to_gltf(pos):
            return [pos[0], pos[2], -pos[1]]
        
        # Map (i, 'ul'|'ur'|'ll'|'lr') -> vertex_index
        vertex_map = {}
        
        num_segments = len(ribbon_vertices)
        for i, cross_section in enumerate(ribbon_vertices):
            # Add lower vertices (always the same)
            ll_pos = to_gltf(cross_section['lower_left'])
            lr_pos = to_gltf(cross_section['lower_right'])
            
            ll_idx = len(positions) // 3
            positions.extend(ll_pos)
            vertex_map[(i, 'll')] = ll_idx
            
            lr_idx = len(positions) // 3
            positions.extend(lr_pos)
            vertex_map[(i, 'lr')] = lr_idx
            
            # Add upper vertices (either multiple with concavity or simple 2-vertex)
            if cross_section.get('has_concavity', False):
                # Concave transversal: multiple vertices across width
                upper_vertices = cross_section['upper_vertices']
                num_upper = len(upper_vertices)
                
                for j, upper_vert in enumerate(upper_vertices):
                    v_pos = to_gltf(upper_vert)
                    v_idx = len(positions) // 3
                    positions.extend(v_pos)
                    
                    # Store mapping
                    if j == 0:
                        vertex_map[(i, 'ul')] = v_idx  # Left edge
                    elif j == num_upper - 1:
                        vertex_map[(i, 'ur')] = v_idx  # Right edge
                    else:
                        vertex_map[(i, f'u{j}')] = v_idx  # Intermediate vertices
                
                # Texture coordinates for upper vertices (vary across width)
                u = i / max(1, num_segments - 1)
                for j in range(num_upper):
                    v_coord = 0.75 + 0.25 * (j / (num_upper - 1))  # 0.75 to 1.0
                    texcoords.extend([u, v_coord])
            else:
                # Flat transversal: simple 2-vertex approach
                ul_pos = to_gltf(cross_section['upper_left'])
                ur_pos = to_gltf(cross_section['upper_right'])
                
                ul_idx = len(positions) // 3
                positions.extend(ul_pos)
                vertex_map[(i, 'ul')] = ul_idx
                
                ur_idx = len(positions) // 3
                positions.extend(ur_pos)
                vertex_map[(i, 'ur')] = ur_idx
                
                # Texture coordinates
                u = i / max(1, num_segments - 1)
                texcoords.extend([u, 0.75])  # upper_left
                texcoords.extend([u, 1.0])   # upper_right
            
            # Texture coordinates for lower vertices (always the same)
            u = i / max(1, num_segments - 1)
            texcoords.extend([u, 0.25])  # lower_left
            texcoords.extend([u, 0.0])   # lower_right
        
        # Create triangles for all surfaces
        for i in range(len(ribbon_vertices) - 1):
            # Get lower vertex indices (always the same)
            ll0 = vertex_map[(i, 'll')]
            lr0 = vertex_map[(i, 'lr')]
            ll1 = vertex_map[(i+1, 'll')]
            lr1 = vertex_map[(i+1, 'lr')]
            
            # Check if either cross-section has concavity
            has_concavity_0 = ribbon_vertices[i].get('has_concavity', False)
            has_concavity_1 = ribbon_vertices[i+1].get('has_concavity', False)
            
            if has_concavity_0 or has_concavity_1:
                # Handle concave transversal: resample cross-sections to same vertex count
                # Get actual 3D positions for both cross-sections
                if has_concavity_0:
                    u0_positions = [ribbon_vertices[i]['upper_vertices'][j] for j in range(ribbon_vertices[i]['num_upper_vertices'])]
                else:
                    u0_positions = [ribbon_vertices[i]['upper_left'], ribbon_vertices[i]['upper_right']]
                
                if has_concavity_1:
                    u1_positions = [ribbon_vertices[i+1]['upper_vertices'][j] for j in range(ribbon_vertices[i+1]['num_upper_vertices'])]
                else:
                    u1_positions = [ribbon_vertices[i+1]['upper_left'], ribbon_vertices[i+1]['upper_right']]
                
                # Resample both to the same count (use maximum to preserve detail)
                # The resampling function already preserves edge vertices exactly
                target_count = max(len(u0_positions), len(u1_positions))
                u0_resampled = self._resample_cross_section_positions(u0_positions, target_count)
                u1_resampled = self._resample_cross_section_positions(u1_positions, target_count)
                
                # Convert to GLTF coordinates and create vertices (or reuse existing ones)
                # Build lookup maps for existing vertices from these cross-sections
                u0_existing_indices = []
                u1_existing_indices = []
                
                if has_concavity_0:
                    for j in range(ribbon_vertices[i]['num_upper_vertices']):
                        if j == 0:
                            u0_existing_indices.append(vertex_map[(i, 'ul')])
                        elif j == ribbon_vertices[i]['num_upper_vertices'] - 1:
                            u0_existing_indices.append(vertex_map[(i, 'ur')])
                        else:
                            u0_existing_indices.append(vertex_map[(i, f'u{j}')])
                else:
                    u0_existing_indices = [vertex_map[(i, 'ul')], vertex_map[(i, 'ur')]]
                
                if has_concavity_1:
                    for j in range(ribbon_vertices[i+1]['num_upper_vertices']):
                        if j == 0:
                            u1_existing_indices.append(vertex_map[(i+1, 'ul')])
                        elif j == ribbon_vertices[i+1]['num_upper_vertices'] - 1:
                            u1_existing_indices.append(vertex_map[(i+1, 'ur')])
                        else:
                            u1_existing_indices.append(vertex_map[(i+1, f'u{j}')])
                else:
                    u1_existing_indices = [vertex_map[(i+1, 'ul')], vertex_map[(i+1, 'ur')]]
                
                # Resample and map to indices (reuse existing vertices when possible)
                # CRITICAL: Always reuse edge vertices exactly to prevent gaps
                u0_indices = []
                u1_indices = []
                tolerance = 1e-6
                
                for j, pos in enumerate(u0_resampled):
                    gltf_pos = to_gltf(pos)
                    
                    # Always reuse edge vertices exactly (first and last)
                    if j == 0 and len(u0_existing_indices) > 0:
                        # First vertex: always use the first existing vertex
                        u0_indices.append(u0_existing_indices[0])
                        continue
                    elif j == len(u0_resampled) - 1 and len(u0_existing_indices) > 0:
                        # Last vertex: always use the last existing vertex
                        u0_indices.append(u0_existing_indices[-1])
                        continue
                    
                    # For intermediate vertices, check if they match existing vertices
                    t = j / (target_count - 1) if target_count > 1 else 0.0
                    source_idx = t * (len(u0_existing_indices) - 1)
                    idx_low = int(source_idx)
                    idx_high = min(idx_low + 1, len(u0_existing_indices) - 1)
                    
                    found = False
                    for check_idx in [idx_low, idx_high]:
                        existing_idx = u0_existing_indices[check_idx]
                        existing_pos = positions[existing_idx*3:(existing_idx+1)*3]
                        if all(abs(gltf_pos[k] - existing_pos[k]) < tolerance for k in range(3)):
                            u0_indices.append(existing_idx)
                            found = True
                            break
                    
                    if not found:
                        # Create new vertex
                        v_idx = len(positions) // 3
                        positions.extend(gltf_pos)
                        u0_indices.append(v_idx)
                        # Add texture coordinate
                        u_coord = i / max(1, num_segments - 1)
                        v_coord = 0.75 + 0.25 * (j / max(1, target_count - 1))
                        texcoords.extend([u_coord, v_coord])
                
                for j, pos in enumerate(u1_resampled):
                    gltf_pos = to_gltf(pos)
                    
                    # Always reuse edge vertices exactly (first and last)
                    if j == 0 and len(u1_existing_indices) > 0:
                        # First vertex: always use the first existing vertex
                        u1_indices.append(u1_existing_indices[0])
                        continue
                    elif j == len(u1_resampled) - 1 and len(u1_existing_indices) > 0:
                        # Last vertex: always use the last existing vertex
                        u1_indices.append(u1_existing_indices[-1])
                        continue
                    
                    # For intermediate vertices, check if they match existing vertices
                    t = j / (target_count - 1) if target_count > 1 else 0.0
                    source_idx = t * (len(u1_existing_indices) - 1)
                    idx_low = int(source_idx)
                    idx_high = min(idx_low + 1, len(u1_existing_indices) - 1)
                    
                    found = False
                    for check_idx in [idx_low, idx_high]:
                        existing_idx = u1_existing_indices[check_idx]
                        existing_pos = positions[existing_idx*3:(existing_idx+1)*3]
                        if all(abs(gltf_pos[k] - existing_pos[k]) < tolerance for k in range(3)):
                            u1_indices.append(existing_idx)
                            found = True
                            break
                    
                    if not found:
                        # Create new vertex
                        v_idx = len(positions) // 3
                        positions.extend(gltf_pos)
                        u1_indices.append(v_idx)
                        # Add texture coordinate
                        u_coord = (i + 1) / max(1, num_segments - 1)
                        v_coord = 0.75 + 0.25 * (j / max(1, target_count - 1))
                        texcoords.extend([u_coord, v_coord])
                
                # Connect resampled vertices 1-to-1 (standard approach for ribbon meshes)
                for j in range(target_count - 1):
                    u0_j = u0_indices[j]
                    u0_j1 = u0_indices[j + 1]
                    u1_j = u1_indices[j]
                    u1_j1 = u1_indices[j + 1]
                    
                    # Top surface: counter-clockwise when viewed from above (normal points up)
                    indices.extend([u0_j, u0_j1, u1_j1])
                    indices.extend([u0_j, u1_j1, u1_j])
                
                # Left side: outward-facing normals (reversed winding)
                ul0 = u0_indices[0]
                ul1 = u1_indices[0]
                indices.extend([ul0, ll1, ll0])
                indices.extend([ul0, ul1, ll1])
                
                # Right side: outward-facing normals (reversed winding)
                ur0 = u0_indices[-1]
                ur1 = u1_indices[-1]
                indices.extend([ur0, lr1, ur1])
                indices.extend([ur0, lr0, lr1])
            else:
                # Both flat: simple 2-vertex connection
                ul0 = vertex_map[(i, 'ul')]
                ur0 = vertex_map[(i, 'ur')]
                ul1 = vertex_map[(i+1, 'ul')]
                ur1 = vertex_map[(i+1, 'ur')]
                
                # Top surface: counter-clockwise when viewed from above
                indices.extend([ul0, ur0, ur1])
                indices.extend([ul0, ur1, ul1])
                
                # Left side: outward-facing normals
                indices.extend([ul0, ll1, ll0])
                indices.extend([ul0, ul1, ll1])
                
                # Right side: outward-facing normals
                indices.extend([ur0, lr1, ur1])
                indices.extend([ur0, lr0, lr1])
            
            # Bottom surface: ordered for visibility from below
            # When viewed from below (positive Y direction), quad should be counter-clockwise: ll0 -> ll1 -> lr1 -> lr0
            indices.extend([ll0, ll1, lr1])  # Triangle 1: counter-clockwise when viewed from below
            indices.extend([ll0, lr1, lr0])  # Triangle 2: counter-clockwise when viewed from below
        
        # Add sloped end caps (similar to original script)
        if len(ribbon_vertices) > 0:
            first_cross = ribbon_vertices[0]
            # Get upper Z - handle both flat and concave cases
            if first_cross.get('has_concavity', False):
                # Use average of left and right edge Z
                upper_verts = first_cross['upper_vertices']
                upper_z = (upper_verts[0][2] + upper_verts[-1][2]) / 2.0
            else:
                upper_z = first_cross['upper_left'][2]
            lower_z = first_cross['lower_left'][2]
            sidewall_height = upper_z - lower_z
            horizontal_offset = sidewall_height
            end_cap_offset = horizontal_offset / math.sqrt(2.0)
            
            def from_gltf(pos):
                return [pos[0], -pos[2], pos[1]]
            
            # Start end cap
            if len(ribbon_vertices) > 1:
                second_cross = ribbon_vertices[1]
                # Calculate center for first cross-section
                if first_cross.get('has_concavity', False):
                    upper_verts_0 = first_cross['upper_vertices']
                    center_first = [
                        (upper_verts_0[0][0] + upper_verts_0[-1][0] + 
                         first_cross['lower_left'][0] + first_cross['lower_right'][0]) / 4.0,
                        (upper_verts_0[0][1] + upper_verts_0[-1][1] + 
                         first_cross['lower_left'][1] + first_cross['lower_right'][1]) / 4.0,
                        (upper_verts_0[0][2] + upper_verts_0[-1][2] + 
                         first_cross['lower_left'][2] + first_cross['lower_right'][2]) / 4.0
                    ]
                else:
                    center_first = [
                        (first_cross['upper_left'][0] + first_cross['upper_right'][0] + 
                         first_cross['lower_left'][0] + first_cross['lower_right'][0]) / 4.0,
                        (first_cross['upper_left'][1] + first_cross['upper_right'][1] + 
                         first_cross['lower_left'][1] + first_cross['lower_right'][1]) / 4.0,
                        (first_cross['upper_left'][2] + first_cross['upper_right'][2] + 
                         first_cross['lower_left'][2] + first_cross['lower_right'][2]) / 4.0
                    ]
                # Calculate center for second cross-section
                if second_cross.get('has_concavity', False):
                    upper_verts_1 = second_cross['upper_vertices']
                    center_second = [
                        (upper_verts_1[0][0] + upper_verts_1[-1][0] + 
                         second_cross['lower_left'][0] + second_cross['lower_right'][0]) / 4.0,
                        (upper_verts_1[0][1] + upper_verts_1[-1][1] + 
                         second_cross['lower_left'][1] + second_cross['lower_right'][1]) / 4.0,
                        (upper_verts_1[0][2] + upper_verts_1[-1][2] + 
                         second_cross['lower_left'][2] + second_cross['lower_right'][2]) / 4.0
                    ]
                else:
                    center_second = [
                        (second_cross['upper_left'][0] + second_cross['upper_right'][0] + 
                         second_cross['lower_left'][0] + second_cross['lower_right'][0]) / 4.0,
                        (second_cross['upper_left'][1] + second_cross['upper_right'][1] + 
                         second_cross['lower_left'][1] + second_cross['lower_right'][1]) / 4.0,
                        (second_cross['upper_left'][2] + second_cross['upper_right'][2] + 
                         second_cross['lower_left'][2] + second_cross['lower_right'][2]) / 4.0
                    ]
                dir_x = center_second[0] - center_first[0]
                dir_y = center_second[1] - center_first[1]
                dir_z = center_second[2] - center_first[2]
                dir_len = math.sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z)
                if dir_len > 1e-9:
                    dir_x /= dir_len
                    dir_y /= dir_len
                    dir_z /= dir_len
                else:
                    dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            else:
                # Fallback: use perpendicular direction
                if first_cross.get('has_concavity', False):
                    upper_verts = first_cross['upper_vertices']
                    dir_x = upper_verts[0][0] - upper_verts[-1][0]
                    dir_y = upper_verts[0][1] - upper_verts[-1][1]
                    dir_z = upper_verts[0][2] - upper_verts[-1][2]
                else:
                    dir_x = first_cross['upper_left'][0] - first_cross['upper_right'][0]
                    dir_y = first_cross['upper_left'][1] - first_cross['upper_right'][1]
                    dir_z = first_cross['upper_left'][2] - first_cross['upper_right'][2]
                dir_len = math.sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z)
                if dir_len > 1e-9:
                    dir_x /= dir_len
                    dir_y /= dir_len
                    dir_z /= dir_len
                else:
                    dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            
            # Reverse direction for start cap
            dir_x = -dir_x
            dir_y = -dir_y
            dir_z = -dir_z
            
            ul0 = vertex_map[(0, 'ul')]
            ur0 = vertex_map[(0, 'ur')]
            ll0 = vertex_map[(0, 'll')]
            lr0 = vertex_map[(0, 'lr')]
            
            ll0_pos = positions[ll0*3:(ll0*3+3)]
            lr0_pos = positions[lr0*3:(lr0*3+3)]
            
            ll0_qgis = from_gltf(ll0_pos)
            lr0_qgis = from_gltf(lr0_pos)
            
            ll0_ext = [ll0_qgis[0] + dir_x * end_cap_offset, 
                       ll0_qgis[1] + dir_y * end_cap_offset, 
                       ll0_qgis[2] + dir_z * end_cap_offset]
            lr0_ext = [lr0_qgis[0] + dir_x * end_cap_offset, 
                       lr0_qgis[1] + dir_y * end_cap_offset, 
                       lr0_qgis[2] + dir_z * end_cap_offset]
            
            ll0_ext_idx = len(positions) // 3
            positions.extend(to_gltf(ll0_ext))
            lr0_ext_idx = len(positions) // 3
            positions.extend(to_gltf(lr0_ext))
            
            texcoords.extend([0.0, 0.25])
            texcoords.extend([0.0, 0.0])
            
            # Start cap triangles
            # Left side face
            indices.extend([ul0, ll0_ext_idx, ll0])
            # Right side face
            indices.extend([ur0, lr0, lr0_ext_idx])
            # Bottom face
            indices.extend([ll0, lr0_ext_idx, lr0])
            indices.extend([ll0, ll0_ext_idx, lr0_ext_idx])
            # End face
            indices.extend([ul0, ur0, lr0_ext_idx])
            indices.extend([ul0, lr0_ext_idx, ll0_ext_idx])
            
            # End cap
            last_idx = len(ribbon_vertices) - 1
            last_cross = ribbon_vertices[last_idx]
            
            if len(ribbon_vertices) > 1:
                prev_cross = ribbon_vertices[last_idx - 1]
                # Calculate center for previous cross-section
                if prev_cross.get('has_concavity', False):
                    upper_verts_prev = prev_cross['upper_vertices']
                    center_prev = [
                        (upper_verts_prev[0][0] + upper_verts_prev[-1][0] + 
                         prev_cross['lower_left'][0] + prev_cross['lower_right'][0]) / 4.0,
                        (upper_verts_prev[0][1] + upper_verts_prev[-1][1] + 
                         prev_cross['lower_left'][1] + prev_cross['lower_right'][1]) / 4.0,
                        (upper_verts_prev[0][2] + upper_verts_prev[-1][2] + 
                         prev_cross['lower_left'][2] + prev_cross['lower_right'][2]) / 4.0
                    ]
                else:
                    center_prev = [
                        (prev_cross['upper_left'][0] + prev_cross['upper_right'][0] + 
                         prev_cross['lower_left'][0] + prev_cross['lower_right'][0]) / 4.0,
                        (prev_cross['upper_left'][1] + prev_cross['upper_right'][1] + 
                         prev_cross['lower_left'][1] + prev_cross['lower_right'][1]) / 4.0,
                        (prev_cross['upper_left'][2] + prev_cross['upper_right'][2] + 
                         prev_cross['lower_left'][2] + prev_cross['lower_right'][2]) / 4.0
                    ]
                # Calculate center for last cross-section
                if last_cross.get('has_concavity', False):
                    upper_verts_last = last_cross['upper_vertices']
                    center_last = [
                        (upper_verts_last[0][0] + upper_verts_last[-1][0] + 
                         last_cross['lower_left'][0] + last_cross['lower_right'][0]) / 4.0,
                        (upper_verts_last[0][1] + upper_verts_last[-1][1] + 
                         last_cross['lower_left'][1] + last_cross['lower_right'][1]) / 4.0,
                        (upper_verts_last[0][2] + upper_verts_last[-1][2] + 
                         last_cross['lower_left'][2] + last_cross['lower_right'][2]) / 4.0
                    ]
                else:
                    center_last = [
                        (last_cross['upper_left'][0] + last_cross['upper_right'][0] + 
                         last_cross['lower_left'][0] + last_cross['lower_right'][0]) / 4.0,
                        (last_cross['upper_left'][1] + last_cross['upper_right'][1] + 
                         last_cross['lower_left'][1] + last_cross['lower_right'][1]) / 4.0,
                        (last_cross['upper_left'][2] + last_cross['upper_right'][2] + 
                         last_cross['lower_left'][2] + last_cross['lower_right'][2]) / 4.0
                    ]
                dir_x = center_last[0] - center_prev[0]
                dir_y = center_last[1] - center_prev[1]
                dir_z = center_last[2] - center_prev[2]
                dir_len = math.sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z)
                if dir_len > 1e-9:
                    dir_x /= dir_len
                    dir_y /= dir_len
                    dir_z /= dir_len
                else:
                    dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            else:
                # Fallback: use perpendicular direction
                if last_cross.get('has_concavity', False):
                    upper_verts = last_cross['upper_vertices']
                    dir_x = upper_verts[0][0] - upper_verts[-1][0]
                    dir_y = upper_verts[0][1] - upper_verts[-1][1]
                    dir_z = upper_verts[0][2] - upper_verts[-1][2]
                else:
                    dir_x = last_cross['upper_left'][0] - last_cross['upper_right'][0]
                    dir_y = last_cross['upper_left'][1] - last_cross['upper_right'][1]
                    dir_z = last_cross['upper_left'][2] - last_cross['upper_right'][2]
                dir_len = math.sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z)
                if dir_len > 1e-9:
                    dir_x /= dir_len
                    dir_y /= dir_len
                    dir_z /= dir_len
                else:
                    dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            
            ul1 = vertex_map[(last_idx, 'ul')]
            ur1 = vertex_map[(last_idx, 'ur')]
            ll1 = vertex_map[(last_idx, 'll')]
            lr1 = vertex_map[(last_idx, 'lr')]
            
            ll1_pos = positions[ll1*3:(ll1*3+3)]
            lr1_pos = positions[lr1*3:(lr1*3+3)]
            
            ll1_qgis = from_gltf(ll1_pos)
            lr1_qgis = from_gltf(lr1_pos)
            
            ll1_ext = [ll1_qgis[0] + dir_x * end_cap_offset, 
                       ll1_qgis[1] + dir_y * end_cap_offset, 
                       ll1_qgis[2] + dir_z * end_cap_offset]
            lr1_ext = [lr1_qgis[0] + dir_x * end_cap_offset, 
                       lr1_qgis[1] + dir_y * end_cap_offset, 
                       lr1_qgis[2] + dir_z * end_cap_offset]
            
            ll1_ext_idx = len(positions) // 3
            positions.extend(to_gltf(ll1_ext))
            lr1_ext_idx = len(positions) // 3
            positions.extend(to_gltf(lr1_ext))
            
            texcoords.extend([1.0, 0.25])
            texcoords.extend([1.0, 0.0])
            
            # End cap triangles
            # Left side face
            indices.extend([ul1, ll1_ext_idx, ll1])
            # Right side face
            indices.extend([ur1, lr1, lr1_ext_idx])
            # Bottom face
            indices.extend([ll1, ll1_ext_idx, lr1])
            indices.extend([lr1, ll1_ext_idx, lr1_ext_idx])
            # End face
            indices.extend([ul1, ur1, lr1_ext_idx])
            indices.extend([ul1, lr1_ext_idx, ll1_ext_idx])
        
        # Export to GLTF/GLB
        self._export_gltf_mesh(
            feature,
            positions,
            indices,
            texcoords,
            export_format,
            project_home,
            texture_file,
            input_layer,
            original_geometry,
            dem_layer,
            ribbon_width,
            z_raise,
            num_vertices,
            feedback
        )
    
    def export_ribbon_mesh_simple(self, feature, geometry, base_z_values, smoothed_z_values,
                                  ribbon_width, export_format, project_home, texture_file=None,
                                  input_layer=None, dem_layer=None, z_raise=1.0, feedback=None):
        """
        Export a simple ribbon mesh without curve banking or concavity.
        Used for "only raise trail" mode.
        """
        if geometry is None or geometry.isEmpty():
            return
        
        # Get vertices from geometry - ONLY use XY coordinates, ignore any Z values
        # All Z values come from DEM sampling (base_z_values, smoothed_z_values)
        vertices = list(geometry.vertices())
        if len(vertices) < 2:
            return
        
        if len(base_z_values) != len(vertices) or len(smoothed_z_values) != len(vertices):
            if feedback:
                feedback.pushWarning(f"Feature {feature.id()}: Mismatch in vertex counts, skipping ribbon export.")
            return
        
        ribbon_vertices = []
        half_width = ribbon_width / 2.0
        bottom_offset = 1.0
        
        # Pre-calculate perpendicular vectors
        segment_perpendiculars = []
        for i in range(len(vertices) - 1):
            v0 = vertices[i]
            v1 = vertices[i + 1]
            dx = v1.x() - v0.x()
            dy = v1.y() - v0.y()
            
            length = math.sqrt(dx*dx + dy*dy)
            if length < 1e-9:
                if len(segment_perpendiculars) > 0:
                    segment_perpendiculars.append(segment_perpendiculars[-1])
                else:
                    segment_perpendiculars.append([1.0, 0.0])
                continue
            
            dx /= length
            dy /= length
            nx = -dy
            ny = dx
            segment_perpendiculars.append([nx, ny])
        
        vertex_perpendiculars = []
        for i in range(len(vertices)):
            if i == 0:
                if len(segment_perpendiculars) > 0:
                    vertex_perpendiculars.append(segment_perpendiculars[0])
                else:
                    vertex_perpendiculars.append([1.0, 0.0])
            elif i == len(vertices) - 1:
                if len(segment_perpendiculars) > 0:
                    vertex_perpendiculars.append(segment_perpendiculars[-1])
                else:
                    vertex_perpendiculars.append([1.0, 0.0])
            else:
                seg_before_idx = i - 1
                seg_after_idx = i
                
                if seg_before_idx < len(segment_perpendiculars) and seg_after_idx < len(segment_perpendiculars):
                    nx1, ny1 = segment_perpendiculars[seg_before_idx]
                    nx2, ny2 = segment_perpendiculars[seg_after_idx]
                    nx_avg = (nx1 + nx2) / 2.0
                    ny_avg = (ny1 + ny2) / 2.0
                    length_avg = math.sqrt(nx_avg*nx_avg + ny_avg*ny_avg)
                    if length_avg > 1e-9:
                        nx_avg /= length_avg
                        ny_avg /= length_avg
                    vertex_perpendiculars.append([nx_avg, ny_avg])
                elif seg_before_idx < len(segment_perpendiculars):
                    vertex_perpendiculars.append(segment_perpendiculars[seg_before_idx])
                elif seg_after_idx < len(segment_perpendiculars):
                    vertex_perpendiculars.append(segment_perpendiculars[seg_after_idx])
                else:
                    vertex_perpendiculars.append([1.0, 0.0])
        
        # Generate vertices
        for i in range(len(vertices)):
            v = vertices[i]
            x = v.x()  # Only use XY from geometry
            y = v.y()  # Only use XY from geometry
            base_z = smoothed_z_values[i]  # Always use DEM-sampled Z values
            # Calculate upper Z from DEM-sampled base elevation + raise
            # Completely ignore any Z values from input linestring geometry
            upper_z = base_z + z_raise
            lower_z = base_z - bottom_offset
            sidewall_height = upper_z - lower_z
            horizontal_offset = sidewall_height
            
            nx, ny = vertex_perpendiculars[i]
            
            upper_left = [x + nx * half_width, y + ny * half_width, upper_z]
            upper_right = [x - nx * half_width, y - ny * half_width, upper_z]
            lower_left = [x + nx * (half_width + horizontal_offset), y + ny * (half_width + horizontal_offset), lower_z]
            lower_right = [x - nx * (half_width + horizontal_offset), y - ny * (half_width + horizontal_offset), lower_z]
            
            ribbon_vertices.append({
                'lower_left': lower_left,
                'lower_right': lower_right,
                'upper_left': upper_left,
                'upper_right': upper_right
            })
        
        if len(ribbon_vertices) < 2:
            return
        
        # Build mesh
        positions = []
        indices = []
        texcoords = []
        
        def to_gltf(pos):
            return [pos[0], pos[2], -pos[1]]
        
        vertex_map = {}
        num_segments = len(ribbon_vertices)
        
        for i, cross_section in enumerate(ribbon_vertices):
            ul_pos = to_gltf(cross_section['upper_left'])
            ur_pos = to_gltf(cross_section['upper_right'])
            ll_pos = to_gltf(cross_section['lower_left'])
            lr_pos = to_gltf(cross_section['lower_right'])
            
            ul_idx = len(positions) // 3
            positions.extend(ul_pos)
            vertex_map[(i, 'ul')] = ul_idx
            
            ur_idx = len(positions) // 3
            positions.extend(ur_pos)
            vertex_map[(i, 'ur')] = ur_idx
            
            ll_idx = len(positions) // 3
            positions.extend(ll_pos)
            vertex_map[(i, 'll')] = ll_idx
            
            lr_idx = len(positions) // 3
            positions.extend(lr_pos)
            vertex_map[(i, 'lr')] = lr_idx
            
            u = i / max(1, num_segments - 1)
            texcoords.extend([u, 0.75])  # upper_left
            texcoords.extend([u, 1.0])   # upper_right
            texcoords.extend([u, 0.25])  # lower_left
            texcoords.extend([u, 0.0])   # lower_right
        
        # Create triangles
        for i in range(len(ribbon_vertices) - 1):
            ul0 = vertex_map[(i, 'ul')]
            ur0 = vertex_map[(i, 'ur')]
            ll0 = vertex_map[(i, 'll')]
            lr0 = vertex_map[(i, 'lr')]
            ul1 = vertex_map[(i+1, 'ul')]
            ur1 = vertex_map[(i+1, 'ur')]
            ll1 = vertex_map[(i+1, 'll')]
            lr1 = vertex_map[(i+1, 'lr')]
            
            # Top surface
            indices.extend([ul0, ur0, ur1])
            indices.extend([ul0, ur1, ul1])
            
            # Bottom surface: ordered for visibility from below
            # When viewed from below (positive Y direction), quad should be counter-clockwise: ll0 -> ll1 -> lr1 -> lr0
            indices.extend([ll0, ll1, lr1])  # Triangle 1: counter-clockwise when viewed from below
            indices.extend([ll0, lr1, lr0])  # Triangle 2: counter-clockwise when viewed from below
            
            # Left side
            indices.extend([ul0, ll1, ll0])
            indices.extend([ul0, ul1, ll1])
            
            # Right side
            indices.extend([ur0, lr1, ur1])
            indices.extend([ur0, lr0, lr1])
        
        # Add end caps (simplified)
        if len(ribbon_vertices) > 0:
            first_cross = ribbon_vertices[0]
            upper_z = first_cross['upper_left'][2]
            lower_z = first_cross['lower_left'][2]
            sidewall_height = upper_z - lower_z
            horizontal_offset = sidewall_height
            end_cap_offset = horizontal_offset / math.sqrt(2.0)
            
            def from_gltf(pos):
                return [pos[0], -pos[2], pos[1]]
            
            # Start cap
            if len(ribbon_vertices) > 1:
                second_cross = ribbon_vertices[1]
                center_first = [
                    (first_cross['upper_left'][0] + first_cross['upper_right'][0] + 
                     first_cross['lower_left'][0] + first_cross['lower_right'][0]) / 4.0,
                    (first_cross['upper_left'][1] + first_cross['upper_right'][1] + 
                     first_cross['lower_left'][1] + first_cross['lower_right'][1]) / 4.0,
                    (first_cross['upper_left'][2] + first_cross['upper_right'][2] + 
                     first_cross['lower_left'][2] + first_cross['lower_right'][2]) / 4.0
                ]
                center_second = [
                    (second_cross['upper_left'][0] + second_cross['upper_right'][0] + 
                     second_cross['lower_left'][0] + second_cross['lower_right'][0]) / 4.0,
                    (second_cross['upper_left'][1] + second_cross['upper_right'][1] + 
                     second_cross['lower_left'][1] + second_cross['lower_right'][1]) / 4.0,
                    (second_cross['upper_left'][2] + second_cross['upper_right'][2] + 
                     second_cross['lower_left'][2] + second_cross['lower_right'][2]) / 4.0
                ]
                dir_x = center_second[0] - center_first[0]
                dir_y = center_second[1] - center_first[1]
                dir_z = center_second[2] - center_first[2]
                dir_len = math.sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z)
                if dir_len > 1e-9:
                    dir_x /= dir_len
                    dir_y /= dir_len
                    dir_z /= dir_len
                else:
                    dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            else:
                dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            
            dir_x = -dir_x
            dir_y = -dir_y
            dir_z = -dir_z
            
            ul0 = vertex_map[(0, 'ul')]
            ur0 = vertex_map[(0, 'ur')]
            ll0 = vertex_map[(0, 'll')]
            lr0 = vertex_map[(0, 'lr')]
            
            ll0_pos = positions[ll0*3:(ll0*3+3)]
            lr0_pos = positions[lr0*3:(lr0*3+3)]
            
            ll0_qgis = from_gltf(ll0_pos)
            lr0_qgis = from_gltf(lr0_pos)
            
            ll0_ext = [ll0_qgis[0] + dir_x * end_cap_offset, 
                       ll0_qgis[1] + dir_y * end_cap_offset, 
                       ll0_qgis[2] + dir_z * end_cap_offset]
            lr0_ext = [lr0_qgis[0] + dir_x * end_cap_offset, 
                       lr0_qgis[1] + dir_y * end_cap_offset, 
                       lr0_qgis[2] + dir_z * end_cap_offset]
            
            ll0_ext_idx = len(positions) // 3
            positions.extend(to_gltf(ll0_ext))
            lr0_ext_idx = len(positions) // 3
            positions.extend(to_gltf(lr0_ext))
            
            texcoords.extend([0.0, 0.25])
            texcoords.extend([0.0, 0.0])
            
            indices.extend([ul0, ll0_ext_idx, ll0])
            indices.extend([ur0, lr0, lr0_ext_idx])
            indices.extend([ll0, lr0_ext_idx, lr0])
            indices.extend([ll0, ll0_ext_idx, lr0_ext_idx])
            indices.extend([ul0, ur0, lr0_ext_idx])
            indices.extend([ul0, lr0_ext_idx, ll0_ext_idx])
            
            # End cap
            last_idx = len(ribbon_vertices) - 1
            last_cross = ribbon_vertices[last_idx]
            
            if len(ribbon_vertices) > 1:
                prev_cross = ribbon_vertices[last_idx - 1]
                center_prev = [
                    (prev_cross['upper_left'][0] + prev_cross['upper_right'][0] + 
                     prev_cross['lower_left'][0] + prev_cross['lower_right'][0]) / 4.0,
                    (prev_cross['upper_left'][1] + prev_cross['upper_right'][1] + 
                     prev_cross['lower_left'][1] + prev_cross['lower_right'][1]) / 4.0,
                    (prev_cross['upper_left'][2] + prev_cross['upper_right'][2] + 
                     prev_cross['lower_left'][2] + prev_cross['lower_right'][2]) / 4.0
                ]
                center_last = [
                    (last_cross['upper_left'][0] + last_cross['upper_right'][0] + 
                     last_cross['lower_left'][0] + last_cross['lower_right'][0]) / 4.0,
                    (last_cross['upper_left'][1] + last_cross['upper_right'][1] + 
                     last_cross['lower_left'][1] + last_cross['lower_right'][1]) / 4.0,
                    (last_cross['upper_left'][2] + last_cross['upper_right'][2] + 
                     last_cross['lower_left'][2] + last_cross['lower_right'][2]) / 4.0
                ]
                dir_x = center_last[0] - center_prev[0]
                dir_y = center_last[1] - center_prev[1]
                dir_z = center_last[2] - center_prev[2]
                dir_len = math.sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z)
                if dir_len > 1e-9:
                    dir_x /= dir_len
                    dir_y /= dir_len
                    dir_z /= dir_len
                else:
                    dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            else:
                dir_x, dir_y, dir_z = 1.0, 0.0, 0.0
            
            ul1 = vertex_map[(last_idx, 'ul')]
            ur1 = vertex_map[(last_idx, 'ur')]
            ll1 = vertex_map[(last_idx, 'll')]
            lr1 = vertex_map[(last_idx, 'lr')]
            
            ll1_pos = positions[ll1*3:(ll1*3+3)]
            lr1_pos = positions[lr1*3:(lr1*3+3)]
            
            ll1_qgis = from_gltf(ll1_pos)
            lr1_qgis = from_gltf(lr1_pos)
            
            ll1_ext = [ll1_qgis[0] + dir_x * end_cap_offset, 
                       ll1_qgis[1] + dir_y * end_cap_offset, 
                       ll1_qgis[2] + dir_z * end_cap_offset]
            lr1_ext = [lr1_qgis[0] + dir_x * end_cap_offset, 
                       lr1_qgis[1] + dir_y * end_cap_offset, 
                       lr1_qgis[2] + dir_z * end_cap_offset]
            
            ll1_ext_idx = len(positions) // 3
            positions.extend(to_gltf(ll1_ext))
            lr1_ext_idx = len(positions) // 3
            positions.extend(to_gltf(lr1_ext))
            
            texcoords.extend([1.0, 0.25])
            texcoords.extend([1.0, 0.0])
            
            indices.extend([ul1, ll1_ext_idx, ll1])
            indices.extend([ur1, lr1, lr1_ext_idx])
            indices.extend([ll1, ll1_ext_idx, lr1])
            indices.extend([lr1, ll1_ext_idx, lr1_ext_idx])
            indices.extend([ul1, ur1, lr1_ext_idx])
            indices.extend([ul1, lr1_ext_idx, ll1_ext_idx])
        
        # Export to GLTF/GLB
        self._export_gltf_mesh(
            feature,
            positions,
            indices,
            texcoords,
            export_format,
            project_home,
            texture_file,
            input_layer,
            geometry,
            dem_layer,
            ribbon_width,
            z_raise,
            len(vertices),
            feedback
        )
    
    def _duplicate_edge_vertices(self, positions, indices, texcoords=None, normal_threshold=0.7, feedback=None):
        """
        Duplicate vertices at edges to preserve sharp edges.
        Vertices shared by faces with significantly different normals are duplicated.
        
        Args:
            positions: List of vertex positions [x, y, z, x, y, z, ...]
            indices: List of triangle indices
            texcoords: Optional list of texture coordinates
            normal_threshold: Dot product threshold below which normals are considered different (0.7 = ~45°)
            feedback: Optional QgsProcessingFeedback for progress
            
        Returns:
            (new_positions, new_indices, new_texcoords, vertex_mapping)
            vertex_mapping: dict mapping old_vertex_idx -> list of new_vertex_indices
        """
        num_vertices = len(positions) // 3
        num_triangles = len(indices) // 3
        
        # Step 1: Calculate face normals
        face_normals = []
        for i in range(0, len(indices), 3):
            if i + 2 >= len(indices):
                break
            v0_idx = indices[i]
            v1_idx = indices[i + 1]
            v2_idx = indices[i + 2]
            
            p0 = [positions[v0_idx*3 + j] for j in range(3)]
            p1 = [positions[v1_idx*3 + j] for j in range(3)]
            p2 = [positions[v2_idx*3 + j] for j in range(3)]
            
            edge1 = [p1[j] - p0[j] for j in range(3)]
            edge2 = [p2[j] - p0[j] for j in range(3)]
            
            normal = [
                edge1[1] * edge2[2] - edge1[2] * edge2[1],
                edge1[2] * edge2[0] - edge1[0] * edge2[2],
                edge1[0] * edge2[1] - edge1[1] * edge2[0]
            ]
            
            length = math.sqrt(sum(n*n for n in normal))
            if length > 1e-9:
                normal = [n / length for n in normal]
            else:
                normal = [0.0, 1.0, 0.0]
            
            face_normals.append((i // 3, normal, [v0_idx, v1_idx, v2_idx]))
        
        # Step 2: For each vertex, find all faces that use it and group by similar normals
        vertex_face_groups = {}  # vertex_idx -> list of (face_idx, normal)
        for face_idx, normal, vertex_indices in face_normals:
            for v_idx in vertex_indices:
                if v_idx not in vertex_face_groups:
                    vertex_face_groups[v_idx] = []
                vertex_face_groups[v_idx].append((face_idx, normal))
        
        # Step 3: Group faces by similar normals for each vertex
        vertex_normal_groups = {}  # vertex_idx -> list of normal groups (each group is a list of face indices)
        for v_idx, face_list in vertex_face_groups.items():
            groups = []
            for face_idx, normal in face_list:
                # Try to find a group with similar normal
                found_group = False
                for group in groups:
                    # Check if normal is similar to any normal in this group
                    for _, group_normal in group:
                        dot_product = sum(normal[i] * group_normal[i] for i in range(3))
                        if dot_product >= normal_threshold:
                            group.append((face_idx, normal))
                            found_group = True
                            break
                    if found_group:
                        break
                
                if not found_group:
                    # Create new group
                    groups.append([(face_idx, normal)])
            
            vertex_normal_groups[v_idx] = groups
        
        # Step 4: Duplicate vertices that have multiple normal groups
        new_positions = []
        new_texcoords = [] if texcoords else None
        vertex_mapping = {}  # old_vertex_idx -> list of new_vertex_indices
        new_vertex_count = 0
        
        for old_v_idx in range(num_vertices):
            groups = vertex_normal_groups.get(old_v_idx, [])
            if len(groups) > 1:
                # Vertex needs duplication - create one copy per normal group
                new_indices_for_vertex = []
                for group in groups:
                    new_v_idx = new_vertex_count
                    new_vertex_count += 1
                    new_indices_for_vertex.append(new_v_idx)
                    
                    # Copy position
                    new_positions.extend(positions[old_v_idx*3:(old_v_idx+1)*3])
                    
                    # Copy texture coordinate if available
                    if texcoords and old_v_idx * 2 + 1 < len(texcoords):
                        new_texcoords.extend(texcoords[old_v_idx*2:(old_v_idx+1)*2])
                    elif texcoords:
                        new_texcoords.extend([0.0, 0.0])
                
                vertex_mapping[old_v_idx] = new_indices_for_vertex
            else:
                # No duplication needed - use single vertex
                new_v_idx = new_vertex_count
                new_vertex_count += 1
                new_positions.extend(positions[old_v_idx*3:(old_v_idx+1)*3])
                
                if texcoords and old_v_idx * 2 + 1 < len(texcoords):
                    new_texcoords.extend(texcoords[old_v_idx*2:(old_v_idx+1)*2])
                elif texcoords:
                    new_texcoords.extend([0.0, 0.0])
                
                vertex_mapping[old_v_idx] = [new_v_idx]
        
        # Step 5: Update indices to use correct duplicate vertices
        new_indices = []
        for face_idx, normal, vertex_indices in face_normals:
            # For each vertex in this face, find which duplicate to use
            new_face_vertices = []
            for old_v_idx in vertex_indices:
                groups = vertex_normal_groups[old_v_idx]
                # Find which group this face belongs to
                target_group_idx = None
                for group_idx, group in enumerate(groups):
                    if any(f_idx == face_idx for f_idx, _ in group):
                        target_group_idx = group_idx
                        break
                
                if target_group_idx is not None:
                    new_v_idx = vertex_mapping[old_v_idx][target_group_idx]
                else:
                    # Fallback: use first duplicate
                    new_v_idx = vertex_mapping[old_v_idx][0]
                
                new_face_vertices.append(new_v_idx)
            
            new_indices.extend(new_face_vertices)
        
        if feedback:
            num_duplicated = sum(1 for v in vertex_mapping.values() if len(v) > 1)
            feedback.pushInfo(f"Duplicated {num_duplicated} edge vertices for sharp edges")
        
        return new_positions, new_indices, new_texcoords, vertex_mapping
    
    def _calculate_normals(self, positions, indices, feedback=None):
        """Calculate vertex normals from triangle indices."""
        num_vertices = len(positions) // 3
        normals = [0.0] * len(positions)
        face_normals = []
        
        # Calculate face normals for each triangle
        for i in range(0, len(indices), 3):
            if i + 2 >= len(indices):
                break
            v0_idx = indices[i]
            v1_idx = indices[i + 1]
            v2_idx = indices[i + 2]
            
            # Get vertex positions
            p0 = [positions[v0_idx*3 + j] for j in range(3)]
            p1 = [positions[v1_idx*3 + j] for j in range(3)]
            p2 = [positions[v2_idx*3 + j] for j in range(3)]
            
            # Calculate edge vectors
            edge1 = [p1[j] - p0[j] for j in range(3)]
            edge2 = [p2[j] - p0[j] for j in range(3)]
            
            # Cross product: edge1 × edge2
            normal = [
                edge1[1] * edge2[2] - edge1[2] * edge2[1],
                edge1[2] * edge2[0] - edge1[0] * edge2[2],
                edge1[0] * edge2[1] - edge1[1] * edge2[0]
            ]
            
            # Normalize
            length = math.sqrt(sum(n*n for n in normal))
            if length > 1e-9:
                normal = [n / length for n in normal]
                # Store triangle index for debugging
                face_normals.append((i // 3, v0_idx, v1_idx, v2_idx, normal, p0, p1, p2))
                
                # Add to vertex normals (will average later)
                for vi in [v0_idx, v1_idx, v2_idx]:
                    for j in range(3):
                        normals[vi*3 + j] += normal[j]
            else:
                # Degenerate triangle, use default up vector
                face_normals.append((i // 3, v0_idx, v1_idx, v2_idx, [0.0, 1.0, 0.0], p0, p1, p2))
                for vi in [v0_idx, v1_idx, v2_idx]:
                    normals[vi*3 + 1] += 1.0
        
        # Normalize vertex normals (average of adjacent faces)
        for i in range(0, len(normals), 3):
            length = math.sqrt(sum(normals[i+j]**2 for j in range(3)))
            if length > 1e-9:
                for j in range(3):
                    normals[i+j] /= length
            else:
                # Default up vector if no valid normal
                normals[i+1] = 1.0
        
        # Log normal statistics
        if feedback:
            # Check top surface normals (should point up, positive Y in GLTF)
            up_normals = 0
            down_normals = 0
            side_normals = 0
            
            for i in range(0, len(normals), 3):
                ny = normals[i+1]  # Y component in GLTF (up/down)
                if ny > 0.5:
                    up_normals += 1
                elif ny < -0.5:
                    down_normals += 1
                else:
                    side_normals += 1
            
            feedback.pushInfo(f"Normal statistics: {up_normals} up, {down_normals} down, {side_normals} side")
            
            # Check first few face normals
            for face_data in face_normals[:5]:
                tri_idx, v0, v1, v2, fn, p0, p1, p2 = face_data
                avg_y = (p0[1] + p1[1] + p2[1]) / 3.0
                feedback.pushInfo(f"Face {tri_idx}: vertices [{v0}, {v1}, {v2}], normal [{fn[0]:.3f}, {fn[1]:.3f}, {fn[2]:.3f}], avg_y={avg_y:.2f}")
            
            # Check for potential winding order issues
            # Analyze by normal direction to identify which surfaces might have issues
            strongly_down = 0
            slightly_down = 0
            strongly_up = 0
            slightly_up = 0
            horizontal = 0
            
            # Track top surface faces that point down (these are problematic)
            top_faces_down = []
            bottom_faces_up = []
            
            for face_data in face_normals:
                tri_idx, v0, v1, v2, fn, p0, p1, p2 = face_data
                ny = fn[1]  # Y component in GLTF (up/down)
                avg_y = (p0[1] + p1[1] + p2[1]) / 3.0  # Average Y position (elevation in GLTF)
                
                if ny < -0.7:  # Strongly pointing down
                    strongly_down += 1
                    # If this face is at high elevation (top surface), it's incorrectly wound
                    if avg_y > 0:  # Top surface should be at positive Y
                        top_faces_down.append((tri_idx, v0, v1, v2, avg_y))
                elif ny < -0.3:  # Slightly pointing down
                    slightly_down += 1
                    if avg_y > 0:
                        top_faces_down.append((tri_idx, v0, v1, v2, avg_y))
                elif ny > 0.7:  # Strongly pointing up (top surface)
                    strongly_up += 1
                    # If this face is at low elevation (bottom surface), it's incorrectly wound
                    if avg_y < -1.0:  # Bottom surface should be at lower Y
                        bottom_faces_up.append((tri_idx, v0, v1, v2, avg_y))
                elif ny > 0.3:  # Slightly pointing up
                    slightly_up += 1
                    if avg_y < -1.0:
                        bottom_faces_up.append((tri_idx, v0, v1, v2, avg_y))
                else:  # Horizontal (side faces)
                    horizontal += 1
            
            feedback.pushInfo(f"Face normal analysis: {strongly_up} strongly up, {slightly_up} slightly up, {horizontal} horizontal, {slightly_down} slightly down, {strongly_down} strongly down")
            
            # Report actual problems: top surface faces pointing down, or bottom faces pointing up
            if top_faces_down:
                feedback.pushWarning(f"Found {len(top_faces_down)} top surface faces with downward-pointing normals (winding order issue). First 5: {top_faces_down[:5]}")
            if bottom_faces_up:
                feedback.pushWarning(f"Found {len(bottom_faces_up)} bottom surface faces with upward-pointing normals (winding order issue). First 5: {bottom_faces_up[:5]}")
            
            # The downward faces are likely correct (bottom surface should point down)
            if strongly_down > 0 and not top_faces_down:
                feedback.pushInfo(f"Note: {strongly_down} downward-pointing faces are likely bottom surface (correct - bottom should point down)")
        else:
            # No feedback provided, initialize empty lists
            top_faces_down = []
            bottom_faces_up = []
        
        return normals, top_faces_down, bottom_faces_up
    
    def _export_gltf_mesh(self, feature, positions, indices, texcoords, export_format, project_home, 
                          texture_file=None, input_layer=None, geometry=None, dem_layer=None, ribbon_width=None,
                          z_raise=1.0, num_vertices=None, feedback=None):
        """Export mesh to GLTF or GLB format."""
        # Generate filename: try "name" field (case insensitive), fallback to FID
        filename_base = None
        try:
            field_names = [f.name() for f in feature.fields()]
            name_field = next((n for n in field_names if n.lower() == 'name'), None)
            if name_field is not None:
                val = feature[name_field]
                if val is not None and str(val).strip():
                    name_str = str(val).strip()
                    safe_chars = []
                    for ch in name_str:
                        if ch.isalnum() or ch in ('_', '-', '.'):
                            safe_chars.append(ch)
                        elif ch.isspace():
                            safe_chars.append('_')
                        else:
                            safe_chars.append('_')
                    filename_base = ''.join(safe_chars).strip('_')
        except Exception:
            pass
        
        if not filename_base:
            feature_id = feature.id()
            filename_base = f"ribbon_feature_{feature_id}"
        
        ext = '.glb' if export_format == 'GLB' else '.gltf'
        filename = f"{filename_base}{ext}"
        output_path = project_home / filename
        
        # Collect metadata
        metadata = {}
        
        if input_layer:
            crs = input_layer.crs()
            if crs.isValid():
                metadata['crs'] = {
                    'authid': crs.authid() if crs.authid() else None,
                    'proj4': crs.toProj4() if crs.toProj4() else None
                }
        
        if geometry and not geometry.isEmpty():
            bbox = geometry.boundingBox()
            if bbox and bbox.width() > 0 and bbox.height() > 0:
                metadata['extent'] = {
                    'xmin': bbox.xMinimum(),
                    'ymin': bbox.yMinimum(),
                    'xmax': bbox.xMaximum(),
                    'ymax': bbox.yMaximum(),
                    'center': {
                        'x': bbox.center().x(),
                        'y': bbox.center().y()
                    }
                }
        
        # Duplicate edge vertices to preserve sharp edges at trail surface edges
        # Threshold 0.8 (~37°) duplicates vertices at sharp edges while avoiding unnecessary duplication
        # This preserves sharp edges (90° angles) without duplicating at smooth transitions
        positions, indices, texcoords, vertex_mapping = self._duplicate_edge_vertices(
            positions, indices, texcoords, normal_threshold=0.8, feedback=feedback
        )
        
        # Calculate normals and identify incorrectly wound faces
        normals, top_faces_down, bottom_faces_up = self._calculate_normals(positions, indices, feedback)
        
        # Fix incorrectly wound top surface faces by reversing their vertex order
        if top_faces_down or bottom_faces_up:
            # Create a set of triangle indices that need to be reversed
            triangles_to_reverse = set()
            for tri_idx, v0, v1, v2, avg_y in top_faces_down:
                triangles_to_reverse.add(tri_idx)
            for tri_idx, v0, v1, v2, avg_y in bottom_faces_up:
                triangles_to_reverse.add(tri_idx)
            
            # Reverse the winding order for problematic triangles
            for tri_idx in triangles_to_reverse:
                base_idx = tri_idx * 3
                if base_idx + 2 < len(indices):
                    # Reverse: swap v1 and v2
                    indices[base_idx + 1], indices[base_idx + 2] = indices[base_idx + 2], indices[base_idx + 1]
            
            if feedback:
                feedback.pushInfo(f"Reversed winding order for {len(triangles_to_reverse)} incorrectly wound triangles")
            
            # Recalculate normals after fixing
            normals, _, _ = self._calculate_normals(positions, indices, None)  # Don't log again
        
        # Convert positions to bytes
        positions_bytes = struct.pack(f'<{len(positions)}f', *positions)
        
        # Always use 32-bit unsigned integers for indices (supports up to ~4.3 billion vertices)
        indices_bytes = struct.pack(f'<{len(indices)}I', *indices)
        indices_component_type = 5125  # UNSIGNED_INT
        
        texcoords_bytes = struct.pack(f'<{len(texcoords)}f', *texcoords) if texcoords else b''
        normals_bytes = struct.pack(f'<{len(normals)}f', *normals)
        
        positions_size = len(positions_bytes)
        indices_size = len(indices_bytes)
        texcoords_size = len(texcoords_bytes)
        normals_size = len(normals_bytes)
        
        # Create GLTF structure
        gltf_asset = {
            "version": "2.0",
            "generator": "TRAILISM QGIS Processing"
        }
        
        if metadata:
            gltf_asset["extras"] = metadata
        
        gltf = {
            "asset": gltf_asset,
            "scene": 0,
            "scenes": [{
                "nodes": [0]
            }],
            "nodes": [{
                "mesh": 0
            }],
            "meshes": [{
                "primitives": [{
                    "attributes": {
                        "POSITION": 0,
                        "NORMAL": 3,
                        "TEXCOORD_0": 2
                    },
                    "indices": 1
                }]
            }],
            "accessors": [
                {
                    "bufferView": 0,
                    "componentType": 5126,
                    "count": len(positions) // 3,
                    "type": "VEC3",
                    "min": [min(positions[i::3]) for i in range(3)],
                    "max": [max(positions[i::3]) for i in range(3)]
                },
                {
                    "bufferView": 1,
                    "componentType": indices_component_type,
                    "count": len(indices),
                    "type": "SCALAR"
                },
                {
                    "bufferView": 2,
                    "componentType": 5126,
                    "count": len(texcoords) // 2,
                    "type": "VEC2"
                },
                {
                    "bufferView": 3,
                    "componentType": 5126,
                    "count": len(normals) // 3,
                    "type": "VEC3"
                }
            ],
            "bufferViews": [
                {
                    "buffer": 0,
                    "byteOffset": 0,
                    "byteLength": positions_size,
                    "target": 34962
                },
                {
                    "buffer": 0,
                    "byteOffset": positions_size,
                    "byteLength": indices_size,
                    "target": 34963
                },
                {
                    "buffer": 0,
                    "byteOffset": positions_size + indices_size,
                    "byteLength": texcoords_size,
                    "target": 34962
                },
                {
                    "buffer": 0,
                    "byteOffset": positions_size + indices_size + texcoords_size,
                    "byteLength": normals_size,
                    "target": 34962
                }
            ],
            "buffers": [{
                "byteLength": positions_size + indices_size + texcoords_size + normals_size
            }],
            "materials": [{
                "name": "Earth",
                "doubleSided": True,  # JSON will serialize True as true
                "pbrMetallicRoughness": {
                    "baseColorFactor": [0.5, 0.4, 0.3, 1.0],
                    "metallicFactor": 0.0,
                    "roughnessFactor": 0.8
                }
            }]
        }
        
        # Handle texture if provided
        if texture_file and Path(texture_file).exists():
            try:
                with open(texture_file, 'rb') as f:
                    texture_data = f.read()
                ext = Path(texture_file).suffix.lower()
                if ext == '.png':
                    texture_mime = 'image/png'
                elif ext in ['.jpg', '.jpeg']:
                    texture_mime = 'image/jpeg'
                else:
                    texture_mime = 'image/png'
                
                texture_idx = len(gltf.get("textures", []))
                image_idx = len(gltf.get("images", []))
                
                if "images" not in gltf:
                    gltf["images"] = []
                if "textures" not in gltf:
                    gltf["textures"] = []
                
                if export_format == 'GLB':
                    texture_base64 = base64.b64encode(texture_data).decode('utf-8')
                    gltf["images"].append({
                        "mimeType": texture_mime,
                        "uri": f"data:{texture_mime};base64,{texture_base64}"
                    })
                else:
                    feature_id = feature.id()
                    texture_filename = f"ribbon_feature_{feature_id}_texture{ext}"
                    texture_path = project_home / texture_filename
                    with open(texture_path, 'wb') as f:
                        f.write(texture_data)
                    gltf["images"].append({
                        "uri": texture_filename
                    })
                
                gltf["textures"].append({
                    "source": image_idx
                })
                
                gltf["materials"][0]["pbrMetallicRoughness"]["baseColorTexture"] = {
                    "index": texture_idx
                }
                
            except Exception as e:
                if feedback:
                    feedback.pushWarning(f"Failed to load texture file: {e}")
        
        gltf["meshes"][0]["primitives"][0]["material"] = 0
        
        if export_format == 'GLB':
            self._write_glb(output_path, gltf, positions_bytes + indices_bytes + texcoords_bytes + normals_bytes)
        else:
            feature_id = feature.id()
            bin_filename = f"ribbon_feature_{feature_id}.bin"
            bin_path = project_home / bin_filename
            gltf["buffers"][0]["uri"] = bin_filename
            
            with open(bin_path, 'wb') as f:
                f.write(positions_bytes + indices_bytes + texcoords_bytes + normals_bytes)
            
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(gltf, f, indent=2)
        
        if feedback:
            feature_id = feature.id()
            feedback.pushInfo(f"Exported ribbon mesh for feature {feature_id} to {output_path}")
    
    def _write_glb(self, output_path, gltf, binary_data):
        """Write GLB file (binary format)."""
        json_str = json.dumps(gltf, separators=(',', ':'))
        json_bytes = json_str.encode('utf-8')
        
        json_padding = (4 - (len(json_bytes) % 4)) % 4
        json_bytes += b' ' * json_padding
        
        bin_padding = (4 - (len(binary_data) % 4)) % 4
        binary_data += b'\x00' * bin_padding
        
        with open(output_path, 'wb') as f:
            f.write(b'glTF')
            f.write(struct.pack('<I', 2))
            f.write(struct.pack('<I', 12 + 8 + len(json_bytes) + 8 + len(binary_data)))
            
            f.write(struct.pack('<I', len(json_bytes)))
            f.write(b'JSON')
            f.write(json_bytes)
            
            f.write(struct.pack('<I', len(binary_data)))
            f.write(b'BIN\x00')
            f.write(binary_data)

