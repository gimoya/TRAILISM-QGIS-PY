## TRAILISM: Trail Shaping

TRAILISM ALGORITHM:\n\nIMPORTANT: All numeric values use period (.) as decimal separator for consistency.\nComma (,) is used only to separate list items.\nNote: The GUI may display default values with comma (,) due to system locale, but always enter values using period (.) as decimal separator.\n\nOverview:\nBuilds a 3D ribbon mesh from a trail centerline: drapes Z from the DEM, optionally raises the outer edge in curves (banking), and can add shape bumps (arc in-run, sigmoid out-run) along the trail. Exports to GLTF/GLB. The input layer is not modified.\n\nCurve detection and banking:\n- The line is densified to 0.5 m spacing along the axis (all vertices at 0.5 m).\n- Curves are detected by segment-based logic: contiguous same-direction turning vertices. A segment is banked only if (1) effective radius ≤ R_max and (2) arc length along the centerline (first to last turning vertex) ≥ min arc length.\n- For each detected curve, banking rises from the first curve vertex (sigmoid on-curve), plateaus through the bend, then falls over 1 m past the last curve vertex.\n- Transversal concavity (optional) adds a shallow concave cross-section in banked sections to mimic a tread shape.\n\nParameters:\n- DEM layer: Required raster for sampling base elevation (Z).\n- Trail raise height: Height above terrain (m). Default: 0.5, minimum: 0.1. Use period (.) for decimals.\n- Selected features only: If checked, process only selected line features.\n- Outer edge raise: Banking height on the outer edge of curves (m). Default: 1.0, minimum: 0.0.\n- Max curve radius (R_max): Bank only bends with effective radius ≤ R_max (m). Larger R_max = gentler bends qualify.\n- Min curve arc length: Bank only if the turning arc along the centerline is at least this long (m). Filters short kinks.\n- Transversal concavity: Depth of cross-sectional concavity in curves (m). Default: 0.3, range 0–1. Use period (.) for decimals.\n- Back slope angle: Angle of the ribbon sidewall/back from horizontal (degrees). Default: 45. Higher = steeper, lower = flatter. Range: 15–90.\n- Tread edge bevel: Chamfer size (m) for the left and right tread edges. 0 = no bevel. Use period (.) for decimal.\n- Ribbon width: Width of exported mesh (m). Default: 2.0, minimum: 0.5. Optional attribute field overrides per feature.\n- Texture file: Optional image (PNG/JPG) for the ribbon surface.\n- Shape distances: Spacing increments between shapes (cumulative along trail). Comma-separated metres; e.g. 8,8,8 places shapes at 8, 16, 24 m. Keywords: "0" = start (first item only), "end" = line end (last item only). Not replicated.\n- Shape heights / table lengths / sigmoid lengths / arc factors: One value replicates to every shape; otherwise one per distance. Empty uses defaults: height 0.5 m, table 0 m, sigmoid 2.5 m, in-run arc 0, landing arc 0. Arc factor: 0 = sigmoid, 1–98 = circular arc, ≥99 = linear ramp. Landing length factor applies when either side is not sigmoid.\n- Export format: GLTF or GLB.\n\nOutput:\n- Ribbon meshes written to the project home directory (or script directory if no project path). Filename from feature "name" field or FID. Meshes are in local coordinates (origin at mesh centroid) for use with Qgis2threejs.\n- One combined GeoJSON (ribbons_points.geojson) with one Point feature per exported ribbon: geometry at mesh centroid; properties x, y, z, rotation, addend for Qgis2threejs; ribbon_path = full path to that feature\

### Algorithm id
`add_trail_shapes`

### Inputs
- **Input line layer** (VectorLayer)
- **Selected features only** (Boolean)
- **DEM layer** (RasterLayer)
- **Trail raise height (m). Range: ≥0.5. Base tread offset above smoothed DEM (before shapes/banking). Fades to 0 at trail ends.** (Number)
- **Outer edge raise (m). Range: ≥0. Banking height on outer berm in curves; 0 = off. Short curves use reduced peak.** (Number)
- **Max curve radius R_max (m). Range: >0. Bank only if effective radius ≤ R_max; gentler bends (larger R) are skipped.** (Number)
- **Min curve arc length (m). Range: ≥0. Bank only if turning arc along centerline (first→last curve vertex) is at least this long.** (Number)
- **Transversal concavity (m). Range: 0–1. Cross-section dish depth in banked curves; 0 = flat tread.** (Number)
- **Back slope angle (°). Range: 15–90. Sidewall angle from horizontal; higher = steeper cut into slope.** (Number)
- **Tread edge bevel (m). Range: ≥0. Chamfer on left/right tread edges; 0 = sharp edge.** (Number)
- **Ribbon width (m). Range: ≥0.5. Trail width when no attribute field is set; field overrides per feature.** (Number)
- **Ribbon width from attribute (optional). Numeric field overrides ribbon width per feature (min 0.5 m).** (Field)
- **Texture image file (optional)** (File)
- **Shape distances (comma-separated, m). Spacing increments — cumulative along trail (8,8,8 → 8, 16, 24). One per shape, not replicated. "0" = start (first only), "end" = line end (last only). Period (.) for decimals.** (String)
- **Shape heights (comma-separated, m). Range: ≥0. Empty = 0.5 m for all. One value replicates; else one per distance.** (String)
- **Shape table lengths (comma-separated, m). Range: ≥0. Empty = 0 m for all. One value replicates; else one per distance.** (String)
- **Sigmoid lengths (comma-separated, m). Range: >0. Empty = 2.5 m for all. Transition length for in-run and landing. One value replicates; else one per distance.** (String)
- **Arc radius factors — in-run (comma-separated). 0 = sigmoid; 1–98 = arc (1 = quarter circle, >1 = gentler); ≥99 = linear ramp. Empty = 0 for all. One value replicates; else one per distance.** (String)
- **Arc radius factors — landing (comma-separated). 0 = sigmoid; 1–98 = arc fall; ≥99 = linear fall. Empty = 0 for all. One value replicates; else one per distance.** (String)
- **Landing length factor. Range: >0. Default: 1.5. Landing transition = in-run length × factor when in-run or landing uses arc; no scale when both are sigmoid (0).** (Number)
- **Ribbon locations GeoJSON output file (optional)** (FileDestination)
- **Export format** (Enum)

### Parameters
- **Input line layer** (VectorLayer)
- **Selected features only** (Boolean, default True)
- **DEM layer** (RasterLayer)
- **Trail raise height (m). Range: ≥0.5. Base tread offset above smoothed DEM (before shapes/banking). Fades to 0 at trail ends.** (Number, default 0.5)
- **Outer edge raise (m). Range: ≥0. Banking height on outer berm in curves; 0 = off. Short curves use reduced peak.** (Number, default 1)
- **Max curve radius R_max (m). Range: >0. Bank only if effective radius ≤ R_max; gentler bends (larger R) are skipped.** (Number, default CURVE_R_MAX_DEFAULT_M)
- **Min curve arc length (m). Range: ≥0. Bank only if turning arc along centerline (first→last curve vertex) is at least this long.** (Number, default CURVE_MIN_ARC_LENGTH_DEFAULT_M)
- **Transversal concavity (m). Range: 0–1. Cross-section dish depth in banked curves; 0 = flat tread.** (Number, default 0.15)
- **Back slope angle (°). Range: 15–90. Sidewall angle from horizontal; higher = steeper cut into slope.** (Number, default 45.0)
- **Tread edge bevel (m). Range: ≥0. Chamfer on left/right tread edges; 0 = sharp edge.** (Number, default 0.02)
- **Ribbon width (m). Range: ≥0.5. Trail width when no attribute field is set; field overrides per feature.** (Number, default 2)
- **Ribbon width from attribute (optional). Numeric field overrides ribbon width per feature (min 0.5 m).** (Field)
- **Texture image file (optional)** (File)
- **Shape distances (comma-separated, m). Spacing increments — cumulative along trail (8,8,8 → 8, 16, 24). One per shape, not replicated. "0" = start (first only), "end" = line end (last only). Period (.) for decimals.** (String, default '')
- **Shape heights (comma-separated, m). Range: ≥0. Empty = 0.5 m for all. One value replicates; else one per distance.** (String, default '')
- **Shape table lengths (comma-separated, m). Range: ≥0. Empty = 0 m for all. One value replicates; else one per distance.** (String, default '')
- **Sigmoid lengths (comma-separated, m). Range: >0. Empty = 2.5 m for all. Transition length for in-run and landing. One value replicates; else one per distance.** (String, default '')
- **Arc radius factors — in-run (comma-separated). 0 = sigmoid; 1–98 = arc (1 = quarter circle, >1 = gentler); ≥99 = linear ramp. Empty = 0 for all. One value replicates; else one per distance.** (String, default '')
- **Arc radius factors — landing (comma-separated). 0 = sigmoid; 1–98 = arc fall; ≥99 = linear fall. Empty = 0 for all. One value replicates; else one per distance.** (String, default '')
- **Landing length factor. Range: >0. Default: 1.5. Landing transition = in-run length × factor when in-run or landing uses arc; no scale when both are sigmoid (0).** (Number, default SHAPE_LANDING_LENGTH_FACTOR_DEFAULT)
- **Ribbon locations GeoJSON output file (optional)** (FileDestination)
- **Export format** (Enum, default 1  # GLB by default)

### Usage
- Python:
```python
import processing
processing.run('add_trail_shapes', {
    'INPUT': ...,
    'SELECTED_FEATURES_ONLY': ...,
    'DEM': ...,
    'Z_RAISE_PARAM': ...,
    'OUTER_EDGE_RAISE': ...,
    'CURVE_R_MAX': ...,
    'CURVE_MIN_ARC_LENGTH': ...,
    'TRANSVERSAL_CONCAVITY': ...,
    'BACK_SLOPE_ANGLE_DEG': ...,
    'TREAD_BEVEL_SIZE': ...,
    'RIBBON_WIDTH': ...,
    'RIBBON_WIDTH_FIELD': ...,
    'TEXTURE_FILE': ...,
    'SHAPE_DISTANCES': ...,
    'SHAPE_HEIGHTS': ...,
    'SHAPE_TABLE_LENGTHS': ...,
    'SIGMOID_LENGTHS': ...,
    'SHAPE_ARC_RADIUS_FACTOR': ...,
    'SHAPE_ARC_RADIUS_FACTOR_LANDING': ...,
    'SHAPE_LANDING_LENGTH_FACTOR': ...,
    'OUTPUT_GEOJSON': ...,
    'EXPORT_FORMAT': ...,
    'OUTPUT': 'memory:'
})
```
