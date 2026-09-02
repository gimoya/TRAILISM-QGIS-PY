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
    QgsProcessingParameterField,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
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

# Curve detection: segment banks only if BOTH (1) effective radius ≤ R_max and
# (2) arc length along centerline (first→last turning vertex) ≥ min arc length.
CURVE_R_MAX_DEFAULT_M = 10
CURVE_MIN_ARC_LENGTH_DEFAULT_M = 5

# Banking envelope: rise on curve from first curve vertex; fall extends past curve end.
# Transition lengths (BANKING_*_LENGTH_M) are measured on the inner ribbon edge; centerline
# span is scaled by R/(R-w) per curve segment so tight turns do not steepen the sigmoid.
BANKING_RISE_LENGTH_M = 2.5
BANKING_FALL_EXTENSION_M = 2
BANKING_SIGMOID_K = 5.0
BANKING_MESH_STEP_M = 0.25   # densify step when banking active (sigmoid needs >2 samples per transition)
BANKING_SHORT_CURVE_LENGTH_M = 6.0
BANKING_SHORT_CURVE_HEIGHT_FACTOR = 0.75

# Tread fade at ends: raise ramps 0->1 over first N vertices and 1->0 over last N; end caps are skipped
TREAD_FADE_VERTICES = 5

# Shape ramps: in-run = circular arc (tangent to level tread); out-run = sigmoid with fixed steepness.
SHAPE_OUTRUN_SIGMOID_K = 4.5  # hardcoded steepness for shape out-run (fall) only
SHAPE_DEFAULT_HEIGHT_M = 0.5
SHAPE_DEFAULT_TABLE_LENGTH_M = 0.0
SHAPE_DEFAULT_SIGMOID_LENGTH_M = 2.5
SHAPE_DEFAULT_ARC_RADIUS_FACTOR = 1.0

# Texture UV: world-space scale (meters per texture tile); 1.0 = one tile per meter along/across
TEXTURE_SCALE_M = 1.0

class AddCurveBankingToLineAlgorithm(QgsProcessingAlgorithm):
    """
    Segment-based curve detection (same-sign turning angles, total angle and curvature threshold).
    For detected curves, raises only the outer edge with a vertex-based sigmoid banking envelope; exports ribbon to GLTF/GLB.
    """

    INPUT = 'INPUT'
    DEM = 'DEM'
    SELECTED_FEATURES_ONLY = 'SELECTED_FEATURES_ONLY'
    Z_RAISE_PARAM = 'Z_RAISE'
    OUTER_EDGE_RAISE = 'OUTER_EDGE_RAISE'
    CURVE_R_MAX = 'CURVE_R_MAX'
    CURVE_MIN_ARC_LENGTH = 'CURVE_MIN_ARC_LENGTH'
    TRANSVERSAL_CONCAVITY = 'TRANSVERSAL_CONCAVITY'
    BACK_SLOPE_ANGLE_DEG = 'BACK_SLOPE_ANGLE_DEG'
    TREAD_BEVEL_SIZE = 'TREAD_BEVEL_SIZE'
    RIBBON_WIDTH = 'RIBBON_WIDTH'
    RIBBON_WIDTH_FIELD = 'RIBBON_WIDTH_FIELD'
    EXPORT_FORMAT = 'EXPORT_FORMAT'
    TEXTURE_FILE = 'TEXTURE_FILE'
    SHAPE_DISTANCES = 'SHAPE_DISTANCES'
    SHAPE_HEIGHTS = 'SHAPE_HEIGHTS'
    SHAPE_TABLE_LENGTHS = 'SHAPE_TABLE_LENGTHS'
    SIGMOID_LENGTHS = 'SIGMOID_LENGTHS'
    SHAPE_ARC_RADIUS_FACTOR = 'SHAPE_ARC_RADIUS_FACTOR'
    OUTPUT_GEOJSON = 'OUTPUT_GEOJSON'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return AddCurveBankingToLineAlgorithm()

    def name(self):
        return 'add_trail_shapes'

    def displayName(self):
        return 'TRAILISM: Trail Shaping '
    
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
            'Overview:\n'
            'Builds a 3D ribbon mesh from a trail centerline: drapes Z from the DEM, optionally raises the outer edge in curves (banking), and can add shape bumps (arc in-run, sigmoid out-run) along the trail. Exports to GLTF/GLB. The input layer is not modified.\n\n'
            'Curve detection and banking:\n'
            '- The line is densified to 0.5 m spacing along the axis (all vertices at 0.5 m).\n'
            '- Curves are detected by segment-based logic: contiguous same-direction turning vertices. A segment is banked only if (1) effective radius ≤ R_max and (2) arc length along the centerline (first to last turning vertex) ≥ min arc length.\n'
            '- For each detected curve, banking rises from the first curve vertex (sigmoid on-curve), plateaus through the bend, then falls over 1 m past the last curve vertex.\n'
            '- Transversal concavity (optional) adds a shallow concave cross-section in banked sections to mimic a tread shape.\n\n'
            'Parameters:\n'
            '- DEM layer: Required raster for sampling base elevation (Z).\n'
            '- Trail raise height: Height above terrain (m). Default: 0.5, minimum: 0.1. Use period (.) for decimals.\n'
            '- Selected features only: If checked, process only selected line features.\n'
            '- Outer edge raise: Banking height on the outer edge of curves (m). Default: 1.0, minimum: 0.0.\n'
            '- Max curve radius (R_max): Bank only bends with effective radius ≤ R_max (m). Larger R_max = gentler bends qualify.\n'
            '- Min curve arc length: Bank only if the turning arc along the centerline is at least this long (m). Filters short kinks.\n'
            '- Transversal concavity: Depth of cross-sectional concavity in curves (m). Default: 0.3, range 0–1. Use period (.) for decimals.\n'
            '- Back slope angle: Angle of the ribbon sidewall/back from horizontal (degrees). Default: 45. Higher = steeper, lower = flatter. Range: 15–90.\n'
            '- Tread edge bevel: Chamfer size (m) for the left and right tread edges. 0 = no bevel. Use period (.) for decimal.\n'
            '- Ribbon width: Width of exported mesh (m). Default: 2.0, minimum: 0.5. Optional attribute field overrides per feature.\n'
            '- Texture file: Optional image (PNG/JPG) for the ribbon surface.\n'
            '- Shape distances: Reference list (one per shape). Comma-separated metres from trail start. Keywords: "0" = start, "end" = line end (last item only). Not replicated.\n'
            '- Shape heights / table lengths / sigmoid lengths / arc factors: One value replicates to every shape; otherwise one per distance. Empty uses defaults: height 0.5 m, table 0 m, sigmoid 2.0 m, arc 1.0. Heights ≥0. Table ≥0 (0 = no plateau). Sigmoid >0. Arc: 0 = sigmoid in-run, ≥1 = circular arc (1 = quarter circle, >1 = gentler).\n'
            '- Export format: GLTF or GLB.\n\n'
            'Output:\n'
            '- Ribbon meshes written to the project home directory (or script directory if no project path). Filename from feature "name" field or FID. Meshes are in local coordinates (origin at mesh centroid) for use with Qgis2threejs.\n'
            '- One combined GeoJSON (ribbons_points.geojson) with one Point feature per exported ribbon: geometry at mesh centroid; properties x, y, z, rotation, addend for Qgis2threejs; ribbon_path = full path to that feature\'s .glb/.gltf; processing_input = JSON string of algorithm inputs for processing scripts. Original layer is not modified.'
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
                'Trail raise height (m). Range: ≥0.5. Base tread offset above smoothed DEM (before shapes/banking). Fades to 0 at trail ends.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
                minValue=0.5
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.OUTER_EDGE_RAISE,
                'Outer edge raise (m). Range: ≥0. Banking height on outer berm in curves; 0 = off. Short curves use reduced peak.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1,
                minValue=0.0,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.CURVE_R_MAX,
                'Max curve radius R_max (m). Range: >0. Bank only if effective radius ≤ R_max; gentler bends (larger R) are skipped.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=CURVE_R_MAX_DEFAULT_M,
                minValue=0.1,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.CURVE_MIN_ARC_LENGTH,
                'Min curve arc length (m). Range: ≥0. Bank only if turning arc along centerline (first→last curve vertex) is at least this long.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=CURVE_MIN_ARC_LENGTH_DEFAULT_M,
                minValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.TRANSVERSAL_CONCAVITY,
                'Transversal concavity (m). Range: 0–1. Cross-section dish depth in banked curves; 0 = flat tread.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.15,
                minValue=0.0,
                maxValue=1.0,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.BACK_SLOPE_ANGLE_DEG,
                'Back slope angle (°). Range: 15–90. Sidewall angle from horizontal; higher = steeper cut into slope.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=45.0,
                minValue=15.0,
                maxValue=90.0,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.TREAD_BEVEL_SIZE,
                'Tread edge bevel (m). Range: ≥0. Chamfer on left/right tread edges; 0 = sharp edge.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.02,
                minValue=0.0,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.RIBBON_WIDTH,
                'Ribbon width (m). Range: ≥0.5. Trail width when no attribute field is set; field overrides per feature.',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=2,
                minValue=0.5,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterField(
                self.RIBBON_WIDTH_FIELD,
                'Ribbon width from attribute (optional). Numeric field overrides ribbon width per feature (min 0.5 m).',
                parentLayerParameterName=self.INPUT,
                type=QgsProcessingParameterField.Numeric,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterFile(
                self.TEXTURE_FILE,
                'Texture image file (optional)',
                behavior=QgsProcessingParameterFile.File,
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_DISTANCES,
                'Shape distances (comma-separated, m). One per shape — reference list, not replicated. Use "0" for start, "end" for line end (last item). Period (.) for decimals.',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_HEIGHTS,
                'Shape heights (comma-separated, m). Range: ≥0. Empty = 0.5 m for all. One value replicates; else one per distance.',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_TABLE_LENGTHS,
                'Shape table lengths (comma-separated, m). Range: ≥0. Empty = 0 m for all. One value replicates; else one per distance.',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SIGMOID_LENGTHS,
                'Sigmoid lengths (comma-separated, m). Range: >0. Empty = 2.0 m for all. Out-run fall; in-run when arc = 0. One value replicates; else one per distance.',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.SHAPE_ARC_RADIUS_FACTOR,
                'Arc radius factors (comma-separated). 0 = sigmoid in-run; ≥1 = arc (1 = quarter circle, >1 = gentler). Empty = 1.0 for all. One value replicates; else one per distance.',
                defaultValue='',
                optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT_GEOJSON,
                'Ribbon locations GeoJSON output file (optional)',
                fileFilter='GeoJSON (*.geojson)',
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

    def _ribbon_width_for_feature(self, feature, ribbon_width_field, default_width):
        """Per-feature ribbon width: from attribute field if set, else default. Min 0.5 m."""
        if ribbon_width_field and feature.fields().indexOf(ribbon_width_field) >= 0:
            val = feature[ribbon_width_field]
            if val is not None and str(val).strip() != '':
                try:
                    return max(0.5, float(val))
                except (TypeError, ValueError):
                    pass
        return max(0.5, float(default_width))

    @staticmethod
    def _parse_comma_list(value_str):
        items = []
        for item in (value_str or '').split(','):
            item = item.strip()
            if item:
                items.append(item)
        return items

    @staticmethod
    def _parse_shape_replicated_values(
        value_str, n_dist, param_name, required=True, empty_default=None, transform=None
    ):
        """
        Parse a comma-separated shape parameter and replicate a single value to match SHAPE_DISTANCES count.
        """
        parsed = AddCurveBankingToLineAlgorithm._parse_comma_list(value_str)
        if not parsed:
            if empty_default is not None:
                return [empty_default] * n_dist
            if required:
                raise ValueError(f"{param_name} cannot be empty when SHAPE_DISTANCES is set.")
            return []
        values = [transform(x) if transform else float(x) for x in parsed]
        if len(values) == 1 and n_dist >= 1:
            values = values * n_dist
        if len(values) != n_dist:
            raise ValueError(
                f"{param_name} has {len(values)} item(s), SHAPE_DISTANCES has {n_dist}. "
                "Provide one value (replicated) or one per distance."
            )
        return values

    @staticmethod
    def _parse_shape_arc_radius_factors(value_str, n_dist):
        return AddCurveBankingToLineAlgorithm._parse_shape_replicated_values(
            value_str,
            n_dist,
            "SHAPE_ARC_RADIUS_FACTOR",
            required=False,
            empty_default=SHAPE_DEFAULT_ARC_RADIUS_FACTOR,
            transform=lambda x: max(0.0, float(x)),
        )

    def _build_processing_input_json(self, parameters, context, input_layer, dem_layer,
                                      ribbon_width_used, ribbon_width_field_name):
        """Build the processing-script input structure (dict) for one exported ribbon.
        Caller should json.dumps(..., indent=2) for the processing_input property.
        """
        export_format_idx = self.parameterAsEnum(parameters, self.EXPORT_FORMAT, context)
        export_format_val = 1 if export_format_idx == 1 else 0  # 1 = GLB, 0 = GLTF
        input_source = input_layer.source() if input_layer else ''
        dem_source = dem_layer.source() if dem_layer else ''
        texture_str = self.parameterAsString(parameters, self.TEXTURE_FILE, context) or None
        ribbon_width_field_val = ribbon_width_field_name if ribbon_width_field_name else None
        shape_distances_str = self.parameterAsString(parameters, self.SHAPE_DISTANCES, context) or ''
        shape_heights_str = self.parameterAsString(parameters, self.SHAPE_HEIGHTS, context) or ''
        shape_table_lengths_str = self.parameterAsString(parameters, self.SHAPE_TABLE_LENGTHS, context) or ''
        sigmoid_lengths_str = self.parameterAsString(parameters, self.SIGMOID_LENGTHS, context) or ''
        shape_arc_radius_factors_str = self.parameterAsString(parameters, self.SHAPE_ARC_RADIUS_FACTOR, context) or ''
        return {
            "area_units": "m2",
            "distance_units": "meters",
            "ellipsoid": "EPSG:7004",
            "inputs": {
                "BACK_SLOPE_ANGLE_DEG": self.parameterAsDouble(parameters, self.BACK_SLOPE_ANGLE_DEG, context) or 45.0,
                "DEM": dem_source,
                "EXPORT_FORMAT": export_format_val,
                "INPUT": input_source,
                "CURVE_R_MAX": self.parameterAsDouble(parameters, self.CURVE_R_MAX, context) or CURVE_R_MAX_DEFAULT_M,
                "CURVE_MIN_ARC_LENGTH": self.parameterAsDouble(parameters, self.CURVE_MIN_ARC_LENGTH, context) or CURVE_MIN_ARC_LENGTH_DEFAULT_M,
                "OUTER_EDGE_RAISE": self.parameterAsDouble(parameters, self.OUTER_EDGE_RAISE, context) or 1.0,
                "RIBBON_WIDTH": ribbon_width_used,
                "RIBBON_WIDTH_FIELD": ribbon_width_field_val,
                "SELECTED_FEATURES_ONLY": self.parameterAsBool(parameters, self.SELECTED_FEATURES_ONLY, context),
                "SHAPE_DISTANCES": shape_distances_str,
                "SHAPE_HEIGHTS": shape_heights_str,
                "SHAPE_TABLE_LENGTHS": shape_table_lengths_str,
                "SIGMOID_LENGTHS": sigmoid_lengths_str,
                "SIGMOID_STEEPNESS_K": SHAPE_OUTRUN_SIGMOID_K,
                "SHAPE_ARC_RADIUS_FACTOR": shape_arc_radius_factors_str,
                "TEXTURE_FILE": texture_str,
                "TRANSVERSAL_CONCAVITY": self.parameterAsDouble(parameters, self.TRANSVERSAL_CONCAVITY, context) or 0.15,
                "TREAD_BEVEL_SIZE": self.parameterAsDouble(parameters, self.TREAD_BEVEL_SIZE, context) or 0.0,
                "Z_RAISE": self.parameterAsDouble(parameters, self.Z_RAISE_PARAM, context) or 0.5,
            }
        }

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
        
        dem_layer = self.parameterAsRasterLayer(parameters, self.DEM, context)
        if dem_layer is None:
            raise QgsProcessingException("DEM layer is required")
        
        # Read basic parameters (always needed)
        z_raise = self.parameterAsDouble(parameters, self.Z_RAISE_PARAM, context)
        z_raise = float(z_raise) if z_raise is not None else 0.5
        if z_raise < 0.1:
            z_raise = 0.1
        
        ribbon_width = self.parameterAsDouble(parameters, self.RIBBON_WIDTH, context)
        ribbon_width = float(ribbon_width) if ribbon_width is not None else 3.0
        ribbon_width_fields = self.parameterAsFields(parameters, self.RIBBON_WIDTH_FIELD, context)
        ribbon_width_field = ribbon_width_fields[0] if ribbon_width_fields else None
        export_format_idx = self.parameterAsEnum(parameters, self.EXPORT_FORMAT, context)
        export_format = 'GLB' if export_format_idx == 1 else 'GLTF'
        texture_file = self.parameterAsString(parameters, self.TEXTURE_FILE, context) or None
        
        # Curve detection: R_max + min arc length along turning vertices
        curve_window_distance = 1.0
        curve_r_max_m = self.parameterAsDouble(parameters, self.CURVE_R_MAX, context)
        curve_r_max_m = float(curve_r_max_m) if curve_r_max_m is not None else CURVE_R_MAX_DEFAULT_M
        if curve_r_max_m <= 0:
            raise QgsProcessingException('Max curve radius R_max must be > 0.')
        curve_min_arc_length_m = self.parameterAsDouble(parameters, self.CURVE_MIN_ARC_LENGTH, context)
        curve_min_arc_length_m = (
            float(curve_min_arc_length_m)
            if curve_min_arc_length_m is not None
            else CURVE_MIN_ARC_LENGTH_DEFAULT_M
        )
        curve_min_arc_length_m = max(0.0, curve_min_arc_length_m)
        outer_edge_raise = self.parameterAsDouble(parameters, self.OUTER_EDGE_RAISE, context)
        outer_edge_raise = float(outer_edge_raise) if outer_edge_raise is not None else 1.0
        transversal_concavity = self.parameterAsDouble(parameters, self.TRANSVERSAL_CONCAVITY, context)
        transversal_concavity = float(transversal_concavity) if transversal_concavity is not None else 0.0
        back_slope_angle_deg = self.parameterAsDouble(parameters, self.BACK_SLOPE_ANGLE_DEG, context)
        back_slope_angle_deg = float(back_slope_angle_deg) if back_slope_angle_deg is not None else 45.0
        back_slope_angle_deg = max(15.0, min(90.0, back_slope_angle_deg))
        if back_slope_angle_deg >= 90.0:
            back_slope_angle_deg = 89.5  # avoid tan(90°) for factor calculation
        
        tread_bevel_size = self.parameterAsDouble(parameters, self.TREAD_BEVEL_SIZE, context)
        tread_bevel_size = max(0.0, float(tread_bevel_size) if tread_bevel_size is not None else 0.0)

        shape_distances_str = self.parameterAsString(parameters, self.SHAPE_DISTANCES, context) or ''
        shape_heights_str = self.parameterAsString(parameters, self.SHAPE_HEIGHTS, context) or ''
        shape_table_lengths_str = self.parameterAsString(parameters, self.SHAPE_TABLE_LENGTHS, context) or ''
        sigmoid_lengths_str = self.parameterAsString(parameters, self.SIGMOID_LENGTHS, context) or ''
        shape_arc_radius_factors_str = self.parameterAsString(parameters, self.SHAPE_ARC_RADIUS_FACTOR, context) or ''
        
        if feedback:
            feedback.pushInfo(f"Trail raise height: {z_raise} meters")
            feedback.pushInfo(
                f"Curve detection: R_max {curve_r_max_m}m, min arc length {curve_min_arc_length_m}m"
            )
            feedback.pushInfo(f"Outer edge raise: {outer_edge_raise} meters")
            feedback.pushInfo(f"Transversal concavity: {transversal_concavity} meters")
            feedback.pushInfo(f"Back slope angle: {back_slope_angle_deg}°")
            if tread_bevel_size > 0:
                feedback.pushInfo(f"Tread edge bevel: {tread_bevel_size} m")

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
        
        # Track which features are exported and collect point features for combined GeoJSON
        exported_features = []
        geojson_features = []

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
                    # Parse shape parameters if provided
                    shape_distances_raw = []
                    shape_heights = []
                    shape_table_lengths = []
                    sigmoid_lengths = []
                    shape_arc_radius_factors = []
                    first_is_zero = False
                    last_is_end = False
                    use_shapes = False
                    
                    def parse_float_list(value_str):
                        return self._parse_comma_list(value_str)

                    if shape_distances_str:
                        try:
                            shape_distances_raw = parse_float_list(shape_distances_str)
                            n_dist = len(shape_distances_raw)
                            if n_dist == 0:
                                msg = f"Feature {feature.id()}: SHAPE_DISTANCES is empty. Provide at least one distance."
                                if feedback:
                                    feedback.pushWarning(msg)
                                raise QgsProcessingException(msg)

                            shape_heights = self._parse_shape_replicated_values(
                                shape_heights_str, n_dist, "SHAPE_HEIGHTS",
                                required=False, empty_default=SHAPE_DEFAULT_HEIGHT_M,
                            )
                            shape_table_lengths = self._parse_shape_replicated_values(
                                shape_table_lengths_str, n_dist, "SHAPE_TABLE_LENGTHS",
                                required=False, empty_default=SHAPE_DEFAULT_TABLE_LENGTH_M,
                            )
                            sigmoid_lengths = self._parse_shape_replicated_values(
                                sigmoid_lengths_str, n_dist, "SIGMOID_LENGTHS",
                                required=False, empty_default=SHAPE_DEFAULT_SIGMOID_LENGTH_M,
                            )
                            shape_arc_radius_factors = self._parse_shape_arc_radius_factors(
                                shape_arc_radius_factors_str, n_dist
                            )

                            use_shapes = True
                            line_length = geom.length()
                            if line_length > 0:
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
                                defaults_used = []
                                if not self._parse_comma_list(shape_heights_str):
                                    defaults_used.append(f"height={SHAPE_DEFAULT_HEIGHT_M}m")
                                if not self._parse_comma_list(shape_table_lengths_str):
                                    defaults_used.append(f"table={SHAPE_DEFAULT_TABLE_LENGTH_M}m")
                                if not self._parse_comma_list(sigmoid_lengths_str):
                                    defaults_used.append(f"sigmoid={SHAPE_DEFAULT_SIGMOID_LENGTH_M}m")
                                if not self._parse_comma_list(shape_arc_radius_factors_str):
                                    defaults_used.append(f"arc={SHAPE_DEFAULT_ARC_RADIUS_FACTOR}")
                                if defaults_used:
                                    feedback.pushInfo(
                                        f"Feature {feature.id()}: shape defaults applied: {', '.join(defaults_used)}"
                                    )
                                feedback.pushInfo(f"Feature {feature.id()}: {len(shape_heights)} shape(s) configured")
                                for i, (h, t, s, r) in enumerate(
                                    zip(shape_heights, shape_table_lengths, sigmoid_lengths, shape_arc_radius_factors)
                                ):
                                    arc_mode = "sigmoid" if r <= 0 else f"arc×{r:.2f}"
                                    feedback.pushInfo(
                                        f"  Shape {i+1}: height={h:.2f}m, table_length={t:.2f}m, "
                                        f"sigmoid_length={s:.2f}m, in-run={arc_mode}"
                                    )
                        except ValueError as e:
                            msg = f"Feature {feature.id()}: {e}"
                            if feedback:
                                feedback.pushWarning(msg)
                            raise QgsProcessingException(msg)
                        except QgsProcessingException:
                            raise
                        except Exception as e:
                            if feedback:
                                feedback.pushWarning(f"Feature {feature.id()}: Error parsing shape parameters - {str(e)}, continuing without shapes")
                            use_shapes = False
                    
                    w = self._ribbon_width_for_feature(feature, ribbon_width_field, ribbon_width)
                    result = self.apply_curve_banking(
                        geom,
                        dem_provider,
                        z_raise,
                        curve_window_distance,
                        curve_r_max_m,
                        curve_min_arc_length_m,
                        outer_edge_raise,
                        w * 0.5,
                        use_shapes,
                        shape_distances_raw if use_shapes else [],
                        shape_heights if use_shapes else [],
                        shape_table_lengths if use_shapes else [],
                        sigmoid_lengths if use_shapes else [],
                        shape_arc_radius_factors if use_shapes else [],
                        first_is_zero,
                        last_is_end,
                        feedback
                    )
                    
                    if isinstance(result, tuple) and len(result) >= 9:
                        original_geom, base_z_values, smoothed_z_values, curve_flags, curve_directions, banking_left, banking_right, distances, shape_z_offsets = result
                    else:
                        feedback.pushWarning(f"Feature {feature.id()}: Failed to process, skipping.")
                        continue
                    
                    if original_geom is None or original_geom.isEmpty():
                        feedback.pushWarning(f"Feature {feature.id()}: Failed to create output geometry, skipping.")
                        continue
                    
                    # Export ribbon mesh
                    if base_z_values is not None and smoothed_z_values is not None:
                        try:
                            result = self.export_ribbon_mesh(
                                feature,
                                original_geom,
                                base_z_values,
                                smoothed_z_values,
                                curve_flags,
                                curve_directions,
                                banking_left,
                                banking_right,
                                distances,
                                shape_z_offsets,
                                w,
                                transversal_concavity,
                                back_slope_angle_deg,
                                tread_bevel_size,
                                export_format,
                                project_home,
                                texture_file,
                                input_layer,
                                dem_layer,
                                z_raise,
                                feedback
                            )
                            if result is not None:
                                ribbon_path, point_feature = result
                                point_feature["properties"]["ribbon_path"] = ribbon_path
                                point_feature["properties"]["processing_input"] = json.dumps(
                                    self._build_processing_input_json(
                                        parameters, context, input_layer, dem_layer, w, ribbon_width_field
                                    ),
                                    indent=2
                                )
                                geojson_features.append(point_feature)
                                exported_features.append(feature.id())
                        except Exception as e:
                            feedback.pushWarning(f"Feature {feature.id()}: Failed to export ribbon mesh - {str(e)}")
                except Exception as e:
                    feedback.pushWarning(f"Feature {feature.id()}: Error processing - {str(e)}, skipping.")
                    continue
                
                feedback.setProgress(int(current * total))
            
            # Write single GeoJSON with all point locations, ribbon_path and processing_input per feature
            if geojson_features:
                combined_geojson = {"type": "FeatureCollection", "features": geojson_features}
                crs = input_layer.crs() if input_layer and input_layer.crs().isValid() else None
                if crs and crs.authid():
                    combined_geojson["crs"] = {
                        "type": "name",
                        "properties": {"name": f"urn:ogc:def:crs:{crs.authid().replace(':', '::')}"},
                    }
                user_geojson_path = self.parameterAsFileOutput(parameters, self.OUTPUT_GEOJSON, context)
                if user_geojson_path:
                    geojson_path = Path(user_geojson_path)
                else:
                    geojson_path = project_home / "ribbons_points.geojson"
                geojson_path.parent.mkdir(parents=True, exist_ok=True)
                with open(geojson_path, 'w', encoding='utf-8') as f:
                    json.dump(combined_geojson, f, indent=2)
                feedback.pushInfo(f"Exported combined point layer for Qgis2threejs to {geojson_path}")

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

    @staticmethod
    def _shape_inrun_arc_height_factor(t, L, R):
        """
        Circular-arc in-run height factor 0->1. Circle tangent to level at s=0, z(s)=R - sqrt(R^2 - s^2).
        Requires R >= L. Returns height factor at fraction t along in-run (s = t*L).
        """
        if L <= 0:
            return 1.0 if t >= 1.0 else 0.0
        R_eff = max(R, L)
        s = min(t * L, R_eff - 1e-9)
        radicand = R_eff * R_eff - s * s
        if radicand <= 0:
            return 1.0 if t >= 1.0 else 0.0
        z_s = R_eff - math.sqrt(radicand)
        denom = R_eff - math.sqrt(max(0.0, R_eff * R_eff - L * L))
        if denom <= 1e-9:
            return 1.0 if t >= 1.0 else (z_s / L if L > 0 else 0.0)
        return max(0.0, min(1.0, z_s / denom))

    def _shape_inrun_height_factor(self, t, L, arc_factor):
        """In-run height factor 0→1: sigmoid when arc_factor ≤ 0, else circular arc."""
        if arc_factor <= 0:
            return self.normalized_sigmoid(t, SHAPE_OUTRUN_SIGMOID_K)
        R = max(arc_factor * L, L)
        return self._shape_inrun_arc_height_factor(t, L, R)

    @staticmethod
    def _segment_total_angle_rad(turning_angles, start_idx, end_idx):
        total_angle_rad = 0.0
        for j in range(start_idx, end_idx + 1):
            if 0 < j < len(turning_angles):
                total_angle_rad += turning_angles[j]
        return total_angle_rad

    @staticmethod
    def _banking_inner_edge_scale(turning_angles, start_idx, end_idx, curve_length, half_width):
        """
        Centerline length factor so BANKING_*_LENGTH_M apply along the inner ribbon edge.
        L_center = L_inner * R / (R - half_width).
        """
        if half_width <= 0 or curve_length <= 1e-9 or not turning_angles:
            return 1.0
        abs_angle = abs(
            AddCurveBankingToLineAlgorithm._segment_total_angle_rad(
                turning_angles, start_idx, end_idx
            )
        )
        if abs_angle < 1e-6:
            return 1.0
        radius_m = curve_length / abs_angle
        radius_m = max(radius_m, half_width + 1e-3)
        return radius_m / (radius_m - half_width)

    def create_curve_banking_envelope(self, curve_flags, distances, outer_edge_raise,
                                      curve_directions=None, turning_angles=None,
                                      half_width=0.0, feedback=None):
        """
        Distance-based banking envelope per curve segment:
          - Rise/fall lengths scaled by R/(R-half_width) so BANKING_*_LENGTH_M apply on the inner edge.
          - Rise: sigmoid 0→peak from first curve vertex over scaled rise length.
          - Plateau: remaining curve vertices at peak.
          - Fall: sigmoid peak→0 over scaled extension past last curve vertex.
          - Left curve (direction > 0) → banking_right; right curve → banking_left.
        Short curves (< 6 m) use 0.5× peak.
        """
        n_vertices = len(curve_flags)
        banking_left = [0.0] * n_vertices
        banking_right = [0.0] * n_vertices
        if curve_directions is None:
            curve_directions = [0.0] * len(curve_flags)
        
        # Build curve segments: split on curve/straight (curve_flags) and on left/right change (curve_directions)
        curve_segments = []
        in_curve = False
        curve_start_idx = None
        
        def _dir_sign(idx):
            d = curve_directions[idx] if idx < len(curve_directions) else 0.0
            return 1 if d > 0 else (-1 if d < 0 else 0)

        for i in range(len(curve_flags)):
            if curve_flags[i] and not in_curve:
                in_curve = True
                curve_start_idx = i
            elif not curve_flags[i] and in_curve:
                in_curve = False
                curve_end_idx = i - 1
                mid_idx = (curve_start_idx + curve_end_idx) // 2
                seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
                curve_segments.append({
                    'start_idx': curve_start_idx,
                    'end_idx': curve_end_idx,
                    'start_dist': distances[curve_start_idx],
                    'end_dist': distances[curve_end_idx],
                    'direction': seg_dir
                })
            elif curve_flags[i] and in_curve and i > curve_start_idx:
                # Back-to-back curves: same curve_flags but direction changed (left <-> right)
                if _dir_sign(curve_start_idx) != _dir_sign(i):
                    curve_end_idx = i - 1
                    mid_idx = (curve_start_idx + curve_end_idx) // 2
                    seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
                    curve_segments.append({
                        'start_idx': curve_start_idx,
                        'end_idx': curve_end_idx,
                        'start_dist': distances[curve_start_idx],
                        'end_dist': distances[curve_end_idx],
                        'direction': seg_dir
                    })
                    curve_start_idx = i
        
        if in_curve:
            curve_end_idx = len(curve_flags) - 1
            mid_idx = (curve_start_idx + curve_end_idx) // 2
            seg_dir = curve_directions[mid_idx] if mid_idx < len(curve_directions) else 0.0
            curve_segments.append({
                'start_idx': curve_start_idx,
                'end_idx': curve_end_idx,
                'start_dist': distances[curve_start_idx],
                'end_dist': distances[curve_end_idx],
                'direction': seg_dir
            })
        
        SHORT_CURVE_LENGTH = BANKING_SHORT_CURVE_LENGTH_M
        SHORT_CURVE_HEIGHT_FACTOR = BANKING_SHORT_CURVE_HEIGHT_FACTOR

        max_inner_edge_scale = 1.0
        for seg in curve_segments:
            start_idx = seg['start_idx']
            end_idx = seg['end_idx']
            if end_idx < start_idx:
                continue
            curve_start_dist = seg['start_dist']
            curve_end_dist = seg['end_dist']
            curve_length = curve_end_dist - curve_start_dist
            # Single-vertex lobe (no ±1 expand in detection): give a nominal span so rise/fall still run
            if curve_length < 1e-9:
                curve_length = BANKING_MESH_STEP_M

            edge_scale = self._banking_inner_edge_scale(
                turning_angles, start_idx, end_idx, curve_length, half_width
            )
            max_inner_edge_scale = max(max_inner_edge_scale, edge_scale)
            rise_length_m = BANKING_RISE_LENGTH_M * edge_scale
            fall_extension_m = BANKING_FALL_EXTENSION_M * edge_scale
            rise_end_dist = curve_start_dist + rise_length_m
            fall_end_dist = curve_end_dist + fall_extension_m

            is_short = curve_length < SHORT_CURVE_LENGTH
            peak_height = (outer_edge_raise * SHORT_CURVE_HEIGHT_FACTOR) if is_short else outer_edge_raise

            if feedback:
                dir_str = "LEFT" if seg['direction'] > 0 else "RIGHT"
                short_note = " (short: 0.5× height)" if is_short else ""
                radius_note = (
                    f", inner-edge scale {edge_scale:.3f} (R/(R-w))"
                    if half_width > 0 else ""
                )
                feedback.pushInfo(
                    f"Banking envelope {dir_str}: curve vertices {start_idx}-{end_idx}, "
                    f"rise {rise_length_m:.2f}m / fall {fall_extension_m:.2f}m on centerline "
                    f"({BANKING_RISE_LENGTH_M:.1f}/{BANKING_FALL_EXTENSION_M:.1f}m on inner edge{radius_note}), "
                    f"peak {peak_height:.2f}m{short_note}"
                )

            for j in range(start_idx, n_vertices):
                d = distances[j]
                if d < curve_start_dist:
                    continue
                if d > fall_end_dist + 1e-9:
                    break

                if d <= curve_end_dist:
                    if d <= rise_end_dist + 1e-9:
                        t = (d - curve_start_dist) / rise_length_m if rise_length_m > 1e-9 else 1.0
                        envelope = peak_height * self.normalized_sigmoid(t, BANKING_SIGMOID_K)
                    else:
                        envelope = peak_height
                else:
                    t = (d - curve_end_dist) / fall_extension_m if fall_extension_m > 1e-9 else 1.0
                    envelope = peak_height * (
                        1.0 - self.normalized_sigmoid(t, BANKING_SIGMOID_K)
                    )

                if seg['direction'] > 0:
                    if envelope > banking_right[j]:
                        banking_right[j] = envelope
                elif seg['direction'] < 0:
                    if envelope > banking_left[j]:
                        banking_left[j] = envelope
        
        return (banking_left, banking_right, max_inner_edge_scale)
    
    def apply_curve_banking(self, geometry, dem_provider, z_raise, 
                           curve_window_distance,
                           curve_r_max_m, curve_min_arc_length_m,
                           outer_edge_raise, ribbon_half_width=0.0,
                           use_shapes=False, shape_distances_raw=None,
                           shape_heights=None, shape_table_lengths=None, sigmoid_lengths=None,
                           shape_arc_radius_factors=None, first_is_zero=False, last_is_end=False, feedback=None):
        """
        Apply curve detection and banking to a linestring.
        Segment-based curve detection (same-sign turning angles, total angle and curvature threshold).
        Banking envelope rise/fall lengths are inner-edge metres, scaled to centerline via R/(R-w).
        """
        # Validate geometry
        if geometry is None or geometry.isEmpty():
            return (geometry, None, None, None, None, None, None, None, [])
        
        # Store original geometry reference - we'll use it directly for XY positions
        original_geometry = geometry
        
        # Densify along axis (finer step when banking) and run segment-based curve detection.
        mesh_step_m = BANKING_MESH_STEP_M if outer_edge_raise > 1e-9 else 0.5
        if feedback and outer_edge_raise > 1e-9:
            feedback.pushInfo(
                f"Banking mesh sampling: {mesh_step_m} m along axis "
                f"(rise {BANKING_RISE_LENGTH_M} m / fall {BANKING_FALL_EXTENSION_M} m on inner edge, "
                f"scaled by R/(R-w) with w={ribbon_half_width:.2f} m)"
            )
        original_geometry, curve_flags, curve_directions, distances, turning_angles = self._densify_adaptive(
            geometry, curve_window_distance,
            curve_r_max_m, curve_min_arc_length_m, mesh_step_m, feedback
        )
        
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
        
        # Step 2: Smooth the base Z values using 20 meter distance-based windows
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
        
        # Curve directions and flags come from _densify_adaptive (no recalculation).
        # angle_sums are no longer needed (peak validation was removed)
        
        # Step 4: Apply shapes if provided
        shape_z_offsets = [0.0] * len(vertices)
        if use_shapes and shape_distances_raw and shape_heights and shape_table_lengths and sigmoid_lengths:
            if not shape_arc_radius_factors:
                shape_arc_radius_factors = [1.0] * len(shape_distances_raw)
            total_length = distances[-1] if distances else 0.0
            shapes = list(zip(
                [float(d) for d in shape_distances_raw],
                shape_heights,
                shape_table_lengths,
                sigmoid_lengths,
                shape_arc_radius_factors,
                range(len(shape_distances_raw)),
            ))
            
            for i, vertex in enumerate(vertices):
                distance_along = distances[i]
                shape_z = 0.0
                
                # Special case: First distance is 0
                if first_is_zero and len(shapes) > 0:
                    first_height = shape_heights[0]
                    first_table_len = shape_table_lengths[0]
                    first_sigmoid_len = sigmoid_lengths[0]
                    first_arc_factor = shape_arc_radius_factors[0]
                    half_table_length = first_table_len / 2.0
                    half_wave_length = half_table_length + first_sigmoid_len
                    if distance_along <= half_wave_length:
                        if distance_along <= half_table_length:
                            shape_z_offset = first_height
                        else:
                            t = (distance_along - half_table_length) / first_sigmoid_len
                            if t <= 1.0:
                                inrun_val = self._shape_inrun_height_factor(
                                    t, first_sigmoid_len, first_arc_factor
                                )
                                shape_z_offset = first_height * (1.0 - inrun_val)
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
                                sigmoid_val = self.normalized_sigmoid(t, SHAPE_OUTRUN_SIGMOID_K)
                                shape_z_offset = last_height * (1.0 - sigmoid_val)
                            else:
                                shape_z_offset = 0.0
                        else:
                            shape_z_offset = 0.0
                        shape_z_offset = max(0.0, shape_z_offset)
                        shape_z += shape_z_offset
                
                # Normal shapes
                for shape_dist, shape_height, table_len, sigmoid_len, arc_factor, orig_idx in shapes:
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
                                inrun_val = self._shape_inrun_height_factor(t, sigmoid_len, arc_factor)
                                shape_z_offset = shape_height * inrun_val
                            else:
                                shape_z_offset = shape_height
                        elif pos_from_wave_start <= sigmoid_len + table_len:
                            shape_z_offset = shape_height
                        else:
                            t = (pos_from_wave_start - sigmoid_len - table_len) / sigmoid_len
                            if t <= 1.0:
                                sigmoid_val = self.normalized_sigmoid(t, SHAPE_OUTRUN_SIGMOID_K)
                                shape_z_offset = shape_height * (1.0 - sigmoid_val)
                            else:
                                shape_z_offset = 0.0
                        
                        shape_z_offset = max(0.0, shape_z_offset)
                        shape_z += shape_z_offset
                
                shape_z_offsets[i] = shape_z
        
        # Create per-side banking envelope (curve-into-curve: both sides can be raised in transition)
        banking_left, banking_right, max_inner_edge_scale = self.create_curve_banking_envelope(
            curve_flags, distances, outer_edge_raise, curve_directions,
            turning_angles, ribbon_half_width, feedback
        )
        
        if feedback:
            n = len(banking_left)
            max_banking = max(max(banking_left or [0]), max(banking_right or [0]))
            num_with_banking = sum(1 for i in range(n) if banking_left[i] > 0 or banking_right[i] > 0)
            num_dual = sum(
                1 for i in range(n)
                if banking_left[i] > 1e-6 and banking_right[i] > 1e-6
            )
            feedback.pushInfo(
                f"Banking heights calculated: max={max_banking:.3f}m, "
                f"{num_with_banking} vertices with banking, {num_dual} dual-side (L and R > 0) out of {n}"
            )
            # Log S-junctions: direction change while still in curve flags
            def _ds(idx):
                d = curve_directions[idx] if idx < len(curve_directions) else 0.0
                return 1 if d > 0 else (-1 if d < 0 else 0)
            for i in range(1, n):
                if not (curve_flags[i] and curve_flags[i - 1]):
                    continue
                if _ds(i - 1) == 0 or _ds(i) == 0 or _ds(i - 1) == _ds(i):
                    continue
                # Sample from junction through fall/rise window
                fall_window_m = BANKING_FALL_EXTENSION_M * max_inner_edge_scale
                j1 = min(n - 1, i + max(1, int(fall_window_m / max(BANKING_MESH_STEP_M, 0.1))))
                dual_in_window = sum(
                    1 for j in range(i, j1 + 1)
                    if banking_left[j] > 1e-6 and banking_right[j] > 1e-6
                )
                feedback.pushInfo(
                    f"S-junction at vertex {i} (dist {distances[i]:.2f}m): "
                    f"dir {_ds(i-1)}→{_ds(i)}, "
                    f"at junction L={banking_left[i]:.3f} R={banking_right[i]:.3f}, "
                    f"dual-side samples in next {fall_window_m:.1f}m: {dual_in_window}"
                )
        
        return (original_geometry, base_z_values, smoothed_z_values, curve_flags,
                curve_directions, banking_left, banking_right, distances, shape_z_offsets)
    
    def _densify_adaptive(self, geometry, curve_window_distance,
                          curve_r_max_m, curve_min_arc_length_m, sample_step_m=0.5, feedback=None):
        """
        Densify geometry to fixed spacing along the axis (default 0.5 m; 0.25 m when banking).
        Process: 1) Densify, 2) Signed turning angles,
        3) Segment-based curve detection: same-sign segments with effective radius ≤ R_max
           and arc length (first→last turning vertex) ≥ min arc length.
        Z values are NOT preserved - they will be sampled from DEM afterwards.
        
        Args:
            geometry: QgsGeometry linestring (2D or 3D, but only X/Y spacing is ensured)
            curve_window_distance: If > 0, curve detection is run.
            curve_r_max_m: Max effective radius (m) for banking.
            curve_min_arc_length_m: Min centerline arc length (m) between first and last turning vertex.
            sample_step_m: Vertex spacing along the line (m).
            feedback: Optional QgsProcessingFeedback for progress
            
        Returns:
            Tuple: (QgsGeometry, curve_flags, curve_directions, distances, turning_angles)
        """
        curve_detection_interval = max(0.1, float(sample_step_m))

        if geometry is None or geometry.isEmpty():
            return (geometry, [], [], [0.0], [])
        
        try:
            geometry_dense = self._densify_at_exact_intervals(geometry, curve_detection_interval, feedback)
            
            # Get vertices at fixed spacing
            vertices_dense = list(geometry_dense.vertices())
            if len(vertices_dense) < 3:
                # Too few vertices, return as-is
                n = len(vertices_dense)
                return (
                    geometry,
                    [False] * n,
                    [0.0] * n,
                    [0.0] * n,
                    [0.0] * n,
                )
            
            # Use exact interval spacing for distances along line
            total_length = geometry.length()
            distances_dense = []
            for i in range(len(vertices_dense)):
                dist = min(i * curve_detection_interval, total_length)
                distances_dense.append(dist)
            
            # Pass 2: Calculate signed turning angles (+180/-180 normalized) for each vertex
            turning_angles = [0.0] * len(vertices_dense)  # Signed angles in radians
            
            for i in range(1, len(vertices_dense) - 1):
                v_prev = vertices_dense[i-1]
                v_curr = vertices_dense[i]
                v_next = vertices_dense[i+1]
                
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
            
            # Same-sign segments: effective radius ≤ R_max and arc length ≥ min (vertex span).
            curve_flags_dense = [False] * len(vertices_dense)
            curve_directions_dense = [0.0] * len(vertices_dense)
            
            if curve_window_distance > 0:
                if feedback:
                    feedback.pushInfo(f"Pass 3: Segment-based curve detection")
                    feedback.pushInfo(
                        f"  R_max: {curve_r_max_m}m, min arc length: {curve_min_arc_length_m}m"
                    )
                
                # Build same-sign segments: contiguous runs of turning_angles with same sign.
                # Bounds are exactly the turning vertices (no ±1 expand into the opposite turn —
                # that overwrote curve_directions at S-junctions and killed dual-side banking).
                segments = []  # list of (start_idx, end_idx, sign)
                run_start = None
                seg_sign = 0
                for i in range(1, len(vertices_dense) - 1):
                    sign = 1 if turning_angles[i] > 0 else (-1 if turning_angles[i] < 0 else 0)
                    if sign != 0:
                        if seg_sign != sign:
                            if run_start is not None:
                                segments.append((run_start, run_end, seg_sign))
                            run_start = i
                            run_end = i
                            seg_sign = sign
                        else:
                            run_end = i
                    else:
                        if run_start is not None:
                            segments.append((run_start, run_end, seg_sign))
                            run_start = None
                        seg_sign = 0
                if run_start is not None:
                    segments.append((run_start, run_end, seg_sign))
                
                curve_count = 0
                n_dist = len(distances_dense)
                for (start_idx, end_idx, direction) in segments:
                    if start_idx > end_idx:
                        continue
                    # Sum turning angles on the same-sign vertices themselves
                    total_angle_rad = self._segment_total_angle_rad(
                        turning_angles, start_idx, end_idx
                    )
                    arc_length = distances_dense[end_idx] - distances_dense[start_idx]
                    if arc_length < 1e-9:
                        continue
                    abs_angle = abs(total_angle_rad)
                    if abs_angle < 1e-6:
                        continue
                    effective_radius_m = arc_length / abs_angle
                    if effective_radius_m <= curve_r_max_m and arc_length >= curve_min_arc_length_m:
                        for j in range(start_idx, end_idx + 1):
                            curve_flags_dense[j] = True
                            curve_directions_dense[j] = float(direction)
                        curve_count += 1
                        if feedback:
                            dir_str = "LEFT" if direction > 0 else "RIGHT"
                            start_dist = distances_dense[start_idx] if start_idx < n_dist else 0.0
                            end_dist = distances_dense[end_idx] if end_idx < n_dist else 0.0
                            feedback.pushInfo(
                                f"  Curve {dir_str}: vertices {start_idx}-{end_idx}, "
                                f"distances {start_dist:.2f}-{end_dist:.2f}m, arc {arc_length:.2f}m, "
                                f"R_eff {effective_radius_m:.2f}m, "
                                f"total angle {math.degrees(abs_angle):.1f}°"
                            )
                
                if feedback:
                    total_curve_vertices = sum(1 for f in curve_flags_dense if f)
                    feedback.pushInfo(f"Pass 3 Summary: Detected {curve_count} curve segments, {total_curve_vertices} curve vertices out of {len(vertices_dense)} total vertices")

            if feedback:
                num_curves = sum(1 for f in curve_flags_dense if f)
                feedback.pushInfo(
                    f"Densification: {len(distances_dense)} vertices at {curve_detection_interval} m "
                    f"along axis (curves: {num_curves} vertices)."
                )
            return (geometry, curve_flags_dense, curve_directions_dense, distances_dense, turning_angles)
                
        except Exception as e:
            if feedback:
                feedback.pushWarning(f"Error in densification: {str(e)}, falling back to uniform {curve_detection_interval}m spacing")
            total_length = geometry.length()
            if total_length > 0:
                fallback_distances = []
                num_intervals = int(math.ceil(total_length / curve_detection_interval))
                for i in range(num_intervals + 1):
                    dist = i * curve_detection_interval
                    if dist <= total_length:
                        fallback_distances.append(dist)
                if fallback_distances[-1] < total_length - 1e-6:
                    fallback_distances.append(total_length)
                n = len(fallback_distances)
                return (
                    geometry,
                    [False] * n,
                    [0.0] * n,
                    fallback_distances,
                    [0.0] * n,
                )
            return (geometry, [], [], [0.0], [])
    
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
                          curve_flags, curve_directions, banking_left, banking_right, distances,
                          shape_z_offsets, ribbon_width, transversal_concavity, back_slope_angle_deg,
                          tread_bevel_size, export_format, project_home, texture_file=None, input_layer=None, dem_layer=None,
                          z_raise=1.0, feedback=None):
        """
        Export a ribbon mesh with curve banking and optional transversal concavity.
        Uses original geometry directly for XY positions to maintain precision.
        For curved sections, raises the outer edge using smooth sigmoid arc envelope.
        Transversal concavity creates a concave (facing down) arc across the width.
        When tread_bevel_size > 0, chamfers the tread edges (flat cross-sections only).
        """
        if original_geometry is None or original_geometry.isEmpty():
            return
        
        if not distances or len(distances) == 0:
            return
        
        num_vertices = len(distances)
        if (len(base_z_values) != num_vertices or len(smoothed_z_values) != num_vertices or
                len(banking_left) != num_vertices or len(banking_right) != num_vertices):
            if feedback:
                feedback.pushWarning(f"Feature {feature.id()}: Mismatch in vertex counts, skipping ribbon export.")
            return
        
        # Generate ribbon mesh vertices
        ribbon_vertices = []
        base_half_width = ribbon_width / 2.0
        bottom_offset = 1.0
        angle_deg = max(15.0, min(89.5, float(back_slope_angle_deg) if back_slope_angle_deg is not None else 45.0))
        back_slope_factor = 1.0 / math.tan(math.radians(angle_deg))
        
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
        
        # Concavity: scalar strength per vertex from max of left/right banking
        ref = max(max(banking_left or [0]), max(banking_right or [0]), 1e-9)
        curve_segment_positions = [min(1.0, max(banking_left[i] if i < len(banking_left) else 0, banking_right[i] if i < len(banking_right) else 0) / ref) for i in range(num_vertices)]
        
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
            
            # Calculate upper Z from DEM-sampled base + raise + shapes (raise faded at first/last N vertices)
            shape_z = shape_z_offsets[i] if i < len(shape_z_offsets) else 0.0
            if TREAD_FADE_VERTICES <= 0:
                fade_factor = 1.0
            else:
                nv = num_vertices
                if i < TREAD_FADE_VERTICES:
                    start_factor = i / max(1, TREAD_FADE_VERTICES - 1)
                else:
                    start_factor = 1.0
                if i >= nv - TREAD_FADE_VERTICES:
                    end_factor = (nv - 1 - i) / max(1, TREAD_FADE_VERTICES - 1)
                else:
                    end_factor = 1.0
                fade_factor = min(start_factor, end_factor)
            effective_raise = z_raise * fade_factor
            upper_z = base_z + effective_raise + shape_z
            
            # Use the perpendicular assigned to this vertex
            nx, ny = vertex_perpendiculars[i]
            
            half_width = base_half_width
            # Constant floor offset; ramp (fade) applies only to tread/raise
            lower_z = base_z - bottom_offset
            
            # Per-side banking: left and right edges raised independently (curve-into-curve: both can be raised)
            bl = banking_left[i] if i < len(banking_left) else 0.0
            br = banking_right[i] if i < len(banking_right) else 0.0
            upper_left_z = upper_z + bl
            upper_right_z = upper_z + br
            
            sidewall_height_left = upper_left_z - lower_z
            sidewall_height_right = upper_right_z - lower_z
            horizontal_offset_left = sidewall_height_left * back_slope_factor
            horizontal_offset_right = sidewall_height_right * back_slope_factor
            lower_left = [x + nx * (half_width + horizontal_offset_left), y + ny * (half_width + horizontal_offset_left), lower_z]
            lower_right = [x - nx * (half_width + horizontal_offset_right), y - ny * (half_width + horizontal_offset_right), lower_z]
            
            # Inner lower edge for concavity clamp: the lower of the two edges when one side is raised
            if bl <= br:
                inner_lower_z = lower_left[2]
            else:
                inner_lower_z = lower_right[2]
            
            if transversal_concavity > 0.0 and curve_segment_positions[i] > 0.0:
                reference_banking = 1.0
                banking_scale = min(1.0, max(bl, br) / reference_banking) if reference_banking > 0 else 0.0
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
                    
                    # Concave sigmoid profile: center is lower (facing down)
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
        scale_m = TEXTURE_SCALE_M

        # Helper function to convert QGIS coords to GLTF coords
        def to_gltf(pos):
            return [pos[0], pos[2], -pos[1]]

        # World-space UV: U = path distance (m), V = lateral position (m); texture repeats every scale_m
        def uv_along(i):
            return distances[i] / scale_m if i < len(distances) else 0.0
        def uv_lateral(lateral_m):
            return lateral_m / scale_m

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
                
                # Texture coordinates for upper vertices (world-space: path distance, lateral m)
                u = uv_along(i)
                for j in range(num_upper):
                    lateral_m = (j / (num_upper - 1)) * ribbon_width if num_upper > 1 else 0.0
                    texcoords.extend([u, uv_lateral(lateral_m)])
                # Bevel at tread edges (curves): same chamfer as flat, using first/last upper vertex
                if tread_bevel_size > 0:
                    ul = upper_vertices[0]
                    ur = upper_vertices[-1]
                    ll = cross_section['lower_left']
                    lr = cross_section['lower_right']
                    eps = 1e-9
                    d_top = [ur[0] - ul[0], ur[1] - ul[1], ur[2] - ul[2]]
                    len_top = math.sqrt(d_top[0]**2 + d_top[1]**2 + d_top[2]**2)
                    ul_bt = [ul[0] + d_top[0] / len_top * tread_bevel_size, ul[1] + d_top[1] / len_top * tread_bevel_size, ul[2] + d_top[2] / len_top * tread_bevel_size] if len_top > eps else list(ul)
                    d_side = [ll[0] - ul[0], ll[1] - ul[1], ll[2] - ul[2]]
                    len_side = math.sqrt(d_side[0]**2 + d_side[1]**2 + d_side[2]**2)
                    ul_bs = [ul[0] + d_side[0] / len_side * tread_bevel_size, ul[1] + d_side[1] / len_side * tread_bevel_size, ul[2] + d_side[2] / len_side * tread_bevel_size] if len_side > eps else list(ul)
                    d_top_r = [ul[0] - ur[0], ul[1] - ur[1], ul[2] - ur[2]]
                    len_top_r = math.sqrt(d_top_r[0]**2 + d_top_r[1]**2 + d_top_r[2]**2)
                    ur_bt = [ur[0] + d_top_r[0] / len_top_r * tread_bevel_size, ur[1] + d_top_r[1] / len_top_r * tread_bevel_size, ur[2] + d_top_r[2] / len_top_r * tread_bevel_size] if len_top_r > eps else list(ur)
                    d_side_r = [lr[0] - ur[0], lr[1] - ur[1], lr[2] - ur[2]]
                    len_side_r = math.sqrt(d_side_r[0]**2 + d_side_r[1]**2 + d_side_r[2]**2)
                    ur_bs = [ur[0] + d_side_r[0] / len_side_r * tread_bevel_size, ur[1] + d_side_r[1] / len_side_r * tread_bevel_size, ur[2] + d_side_r[2] / len_side_r * tread_bevel_size] if len_side_r > eps else list(ur)
                    for pos, key, v_lateral in [(ul_bt, 'ul_bevel_top', 0.1 * ribbon_width), (ul_bs, 'ul_bevel_side', 0.0), (ur_bt, 'ur_bevel_top', 0.9 * ribbon_width), (ur_bs, 'ur_bevel_side', ribbon_width)]:
                        idx = len(positions) // 3
                        positions.extend(to_gltf(pos))
                        vertex_map[(i, key)] = idx
                        texcoords.extend([u, uv_lateral(v_lateral)])
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
                
                # Texture coordinates (world-space)
                u = uv_along(i)
                texcoords.extend([u, uv_lateral(0.0)])       # upper_left
                texcoords.extend([u, uv_lateral(ribbon_width)])  # upper_right

                # Bevel vertices for tread edges (chamfer)
                if tread_bevel_size > 0:
                    ul = cross_section['upper_left']
                    ur = cross_section['upper_right']
                    ll = cross_section['lower_left']
                    lr = cross_section['lower_right']
                    eps = 1e-9
                    # Left: ul_bevel_top (inward along tread), ul_bevel_side (down sidewall)
                    d_top = [ur[0] - ul[0], ur[1] - ul[1], ur[2] - ul[2]]
                    len_top = math.sqrt(d_top[0]**2 + d_top[1]**2 + d_top[2]**2)
                    ul_bt = [ul[0] + d_top[0] / len_top * tread_bevel_size, ul[1] + d_top[1] / len_top * tread_bevel_size, ul[2] + d_top[2] / len_top * tread_bevel_size] if len_top > eps else list(ul)
                    d_side = [ll[0] - ul[0], ll[1] - ul[1], ll[2] - ul[2]]
                    len_side = math.sqrt(d_side[0]**2 + d_side[1]**2 + d_side[2]**2)
                    ul_bs = [ul[0] + d_side[0] / len_side * tread_bevel_size, ul[1] + d_side[1] / len_side * tread_bevel_size, ul[2] + d_side[2] / len_side * tread_bevel_size] if len_side > eps else list(ul)
                    # Right: ur_bevel_top (inward), ur_bevel_side (down)
                    d_top_r = [ul[0] - ur[0], ul[1] - ur[1], ul[2] - ur[2]]
                    len_top_r = math.sqrt(d_top_r[0]**2 + d_top_r[1]**2 + d_top_r[2]**2)
                    ur_bt = [ur[0] + d_top_r[0] / len_top_r * tread_bevel_size, ur[1] + d_top_r[1] / len_top_r * tread_bevel_size, ur[2] + d_top_r[2] / len_top_r * tread_bevel_size] if len_top_r > eps else list(ur)
                    d_side_r = [lr[0] - ur[0], lr[1] - ur[1], lr[2] - ur[2]]
                    len_side_r = math.sqrt(d_side_r[0]**2 + d_side_r[1]**2 + d_side_r[2]**2)
                    ur_bs = [ur[0] + d_side_r[0] / len_side_r * tread_bevel_size, ur[1] + d_side_r[1] / len_side_r * tread_bevel_size, ur[2] + d_side_r[2] / len_side_r * tread_bevel_size] if len_side_r > eps else list(ur)
                    for pos, key, v_lat in [(ul_bt, 'ul_bevel_top', 0.1 * ribbon_width), (ul_bs, 'ul_bevel_side', 0.0), (ur_bt, 'ur_bevel_top', 0.9 * ribbon_width), (ur_bs, 'ur_bevel_side', ribbon_width)]:
                        idx = len(positions) // 3
                        positions.extend(to_gltf(pos))
                        vertex_map[(i, key)] = idx
                        texcoords.extend([u, uv_lateral(v_lat)])
            
            # Texture coordinates for lower vertices (world-space)
            u = uv_along(i)
            texcoords.extend([u, uv_lateral(0.0)])           # lower_left
            texcoords.extend([u, uv_lateral(ribbon_width)])  # lower_right
        
        # Create triangles for all surfaces
        # Cache resampled indices for flat sections so adjacent segments share the same vertices (no sliver at concave/flat transition)
        resampled_flat_cache = {}
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
                # Reuse resampled flat section from previous segment so both segments share the same ring (no sliver)
                if not has_concavity_0 and i in resampled_flat_cache and len(resampled_flat_cache[i]) == target_count:
                    u0_indices = resampled_flat_cache[i]
                else:
                    for j, pos in enumerate(u0_resampled):
                        gltf_pos = to_gltf(pos)
                        if j == 0 and len(u0_existing_indices) > 0:
                            u0_indices.append(u0_existing_indices[0])
                            continue
                        elif j == len(u0_resampled) - 1 and len(u0_existing_indices) > 0:
                            u0_indices.append(u0_existing_indices[-1])
                            continue
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
                            v_idx = len(positions) // 3
                            positions.extend(gltf_pos)
                            u0_indices.append(v_idx)
                            u_coord = uv_along(i)
                            lateral_m = (j / max(1, target_count - 1)) * ribbon_width
                            texcoords.extend([u_coord, uv_lateral(lateral_m)])
                if not has_concavity_0:
                    resampled_flat_cache[i] = u0_indices
                # Reuse resampled flat section from next segment so shared edge uses same vertices
                if not has_concavity_1 and (i + 1) in resampled_flat_cache and len(resampled_flat_cache[i + 1]) == target_count:
                    u1_indices = resampled_flat_cache[i + 1]
                else:
                    for j, pos in enumerate(u1_resampled):
                        gltf_pos = to_gltf(pos)
                        if j == 0 and len(u1_existing_indices) > 0:
                            u1_indices.append(u1_existing_indices[0])
                            continue
                        elif j == len(u1_resampled) - 1 and len(u1_existing_indices) > 0:
                            u1_indices.append(u1_existing_indices[-1])
                            continue
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
                            v_idx = len(positions) // 3
                            positions.extend(gltf_pos)
                            u1_indices.append(v_idx)
                            u_coord = uv_along(i + 1)
                            lateral_m = (j / max(1, target_count - 1)) * ribbon_width
                            texcoords.extend([u_coord, uv_lateral(lateral_m)])
                    if not has_concavity_1:
                        resampled_flat_cache[i + 1] = u1_indices
                # Use bevel at tread edges when tread_bevel_size > 0 (both flat and concave sections have bevel vertices)
                use_concavity_bevel = tread_bevel_size > 0
                if use_concavity_bevel:
                    u0_bevel = [vertex_map[(i, 'ul_bevel_top')]] + u0_indices[1:-1] + [vertex_map[(i, 'ur_bevel_top')]]
                    u1_bevel = [vertex_map[(i+1, 'ul_bevel_top')]] + u1_indices[1:-1] + [vertex_map[(i+1, 'ur_bevel_top')]]
                    u0_top, u1_top = u0_bevel, u1_bevel
                    ul_bs0, ul_bs1 = vertex_map[(i, 'ul_bevel_side')], vertex_map[(i+1, 'ul_bevel_side')]
                    ur_bs0, ur_bs1 = vertex_map[(i, 'ur_bevel_side')], vertex_map[(i+1, 'ur_bevel_side')]
                    ul0, ur0 = u0_indices[0], u0_indices[-1]
                else:
                    u0_top, u1_top = u0_indices, u1_indices
                    ul_bs0, ul_bs1 = u0_indices[0], u1_indices[0]
                    ur_bs0, ur_bs1 = u0_indices[-1], u1_indices[-1]
                    ul0, ur0 = u0_indices[0], u0_indices[-1]
                # Connect resampled vertices 1-to-1 (standard approach for ribbon meshes)
                for j in range(target_count - 1):
                    u0_j = u0_top[j]
                    u0_j1 = u0_top[j + 1]
                    u1_j = u1_top[j]
                    u1_j1 = u1_top[j + 1]
                    # Top surface: counter-clockwise when viewed from above (normal points up)
                    indices.extend([u0_j, u0_j1, u1_j1])
                    indices.extend([u0_j, u1_j1, u1_j])
                # Left side: outward-facing normals (reversed winding)
                indices.extend([ul_bs0, ll1, ll0])
                indices.extend([ul_bs0, ul_bs1, ll1])
                # Right side: outward-facing normals (reversed winding)
                indices.extend([ur_bs0, lr1, ur_bs1])
                indices.extend([ur_bs0, lr0, lr1])
                if use_concavity_bevel:
                    ul_bt0 = vertex_map[(i, 'ul_bevel_top')]
                    ul_bt1 = vertex_map[(i+1, 'ul_bevel_top')]
                    ur_bt0 = vertex_map[(i, 'ur_bevel_top')]
                    ur_bt1 = vertex_map[(i+1, 'ur_bevel_top')]
                    # Bevel as quad between bevel edges only (no corner): top0-top1-side1-side0
                    indices.extend([ul_bt0, ul_bt1, ul_bs1])
                    indices.extend([ul_bt0, ul_bs1, ul_bs0])
                    indices.extend([ur_bt0, ur_bs0, ur_bs1])
                    indices.extend([ur_bt0, ur_bs1, ur_bt1])
            else:
                # Both flat: simple 2-vertex connection (with optional tread bevel)
                ul0 = vertex_map[(i, 'ul')]
                ur0 = vertex_map[(i, 'ur')]
                ul1 = vertex_map[(i+1, 'ul')]
                ur1 = vertex_map[(i+1, 'ur')]
                
                if tread_bevel_size > 0:
                    ul_bt0 = vertex_map[(i, 'ul_bevel_top')]
                    ul_bs0 = vertex_map[(i, 'ul_bevel_side')]
                    ur_bt0 = vertex_map[(i, 'ur_bevel_top')]
                    ur_bs0 = vertex_map[(i, 'ur_bevel_side')]
                    ul_bt1 = vertex_map[(i+1, 'ul_bevel_top')]
                    ul_bs1 = vertex_map[(i+1, 'ul_bevel_side')]
                    ur_bt1 = vertex_map[(i+1, 'ur_bevel_top')]
                    ur_bs1 = vertex_map[(i+1, 'ur_bevel_side')]
                    # Top: boundary is bevel-top vertices
                    indices.extend([ul_bt0, ur_bt0, ur_bt1])
                    indices.extend([ul_bt0, ur_bt1, ul_bt1])
                    # Left side: top of sidewall is bevel-side
                    indices.extend([ul_bs0, ll1, ll0])
                    indices.extend([ul_bs0, ul_bs1, ll1])
                    # Right side
                    indices.extend([ur_bs0, lr1, ur_bs1])
                    indices.extend([ur_bs0, lr0, lr1])
                    # Bevel as quad between bevel edges only (no corner): top0-top1-side1-side0
                    indices.extend([ul_bt0, ul_bt1, ul_bs1])
                    indices.extend([ul_bt0, ul_bs1, ul_bs0])
                    indices.extend([ur_bt0, ur_bs0, ur_bs1])
                    indices.extend([ur_bt0, ur_bs1, ur_bt1])
                else:
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
        
        # Add sloped end caps only when tread fade is off (fade uses first/last N vertices, so no caps)
        if len(ribbon_vertices) > 0 and TREAD_FADE_VERTICES <= 0:
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
            horizontal_offset = sidewall_height * back_slope_factor
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
            if tread_bevel_size > 0 and (0, 'ul_bevel_top') in vertex_map:
                ul0_top = vertex_map[(0, 'ul_bevel_top')]
                ur0_top = vertex_map[(0, 'ur_bevel_top')]
                ul0_side = vertex_map[(0, 'ul_bevel_side')]
                ur0_side = vertex_map[(0, 'ur_bevel_side')]
            else:
                ul0_top, ur0_top = ul0, ur0
                ul0_side, ur0_side = ul0, ur0
            
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

            u_start = uv_along(0)
            texcoords.extend([u_start, uv_lateral(0.0)])
            texcoords.extend([u_start, uv_lateral(ribbon_width)])

            # Start cap triangles (use bevel vertices when tread_bevel_size > 0)
            # Left side face
            indices.extend([ul0_side, ll0_ext_idx, ll0])
            # Right side face
            indices.extend([ur0_side, lr0, lr0_ext_idx])
            # Bottom face
            indices.extend([ll0, lr0_ext_idx, lr0])
            indices.extend([ll0, ll0_ext_idx, lr0_ext_idx])
            # End face (frontal)
            indices.extend([ul0_top, ur0_top, lr0_ext_idx])
            indices.extend([ul0_top, lr0_ext_idx, ll0_ext_idx])
            if tread_bevel_size > 0 and (0, 'ul_bevel_top') in vertex_map:
                # Chamfer fill between end face and side faces
                indices.extend([ul0_top, ul0_side, ll0])
                indices.extend([ur0_top, ur0_side, lr0])
            
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
            if tread_bevel_size > 0 and (last_idx, 'ul_bevel_top') in vertex_map:
                ul1_top = vertex_map[(last_idx, 'ul_bevel_top')]
                ur1_top = vertex_map[(last_idx, 'ur_bevel_top')]
                ul1_side = vertex_map[(last_idx, 'ul_bevel_side')]
                ur1_side = vertex_map[(last_idx, 'ur_bevel_side')]
            else:
                ul1_top, ur1_top = ul1, ur1
                ul1_side, ur1_side = ul1, ur1
            
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

            u_end = uv_along(len(ribbon_vertices) - 1)
            texcoords.extend([u_end, uv_lateral(0.0)])
            texcoords.extend([u_end, uv_lateral(ribbon_width)])

            # End cap triangles (use bevel vertices when tread_bevel_size > 0)
            # Left side face
            indices.extend([ul1_side, ll1_ext_idx, ll1])
            # Right side face
            indices.extend([ur1_side, lr1, lr1_ext_idx])
            # Bottom face
            indices.extend([ll1, ll1_ext_idx, lr1])
            indices.extend([lr1, ll1_ext_idx, lr1_ext_idx])
            # End face (frontal)
            indices.extend([ul1_top, ur1_top, lr1_ext_idx])
            indices.extend([ul1_top, lr1_ext_idx, ll1_ext_idx])
            if tread_bevel_size > 0 and (last_idx, 'ul_bevel_top') in vertex_map:
                # Chamfer fill between end face and side faces
                indices.extend([ul1_top, ul1_side, ll1])
                indices.extend([ur1_top, ur1_side, lr1])
        
        # Export to GLTF/GLB; return (output_path, point_feature) for combined GeoJSON
        return self._export_gltf_mesh(
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

        # Transform to local coordinates (origin at mesh centroid) for Qgis2threejs "3D Model on a point"
        pos_min = [min(positions[i::3]) for i in range(3)]
        pos_max = [max(positions[i::3]) for i in range(3)]
        cx_gltf = (pos_min[0] + pos_max[0]) / 2.0
        cy_gltf = (pos_min[1] + pos_max[1]) / 2.0
        cz_gltf = (pos_min[2] + pos_max[2]) / 2.0
        center = (cx_gltf, cy_gltf, cz_gltf)
        for i in range(len(positions)):
            positions[i] -= center[i % 3]
        half_x = (pos_max[0] - pos_min[0]) / 2.0
        half_y = (pos_max[1] - pos_min[1]) / 2.0
        half_z = (pos_max[2] - pos_min[2]) / 2.0
        position_min_local = [-half_x, -half_y, -half_z]
        position_max_local = [half_x, half_y, half_z]
        # CRS centroid for GeoJSON point (glTF: X=easting, Y=elevation, Z=-northing)
        cx_crs = cx_gltf
        cy_crs = -cz_gltf
        cz_crs = cy_gltf
        if metadata:
            metadata['crs_centroid'] = {'x': cx_crs, 'y': cy_crs, 'z': cz_crs}

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
                    "min": position_min_local,
                    "max": position_max_local
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

        # Point feature for combined GeoJSON (written once after all features by processAlgorithm).
        # Geometry and properties carry full 3D position for Qgis2threejs "3D Model on a point".
        point_feature = {
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [cx_crs, cy_crs, cz_crs]},
            "properties": {
                "x": float(cx_crs),
                "y": float(cy_crs),
                "z": float(cz_crs),
                "rotation": 0.0,
                "addend": 0.0,
            },
        }
        if feedback:
            feature_id = feature.id()
            feedback.pushInfo(f"Exported ribbon mesh for feature {feature_id} to {output_path}")
        return (str(output_path), point_feature)
    
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

