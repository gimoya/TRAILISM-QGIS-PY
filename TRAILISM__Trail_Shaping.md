## TRAILISM: Trail Shaping

### Algorithm id
`add_trail_shapes`

### Inputs
- **Input line layer** (VectorLayer)
- **Selected features only** (Boolean)
- **DEM layer** (RasterLayer)
- **Trail raise height (meters). Default: 0.5, minimum: 0.1. Use period (.) for decimal separator** (Number)
- **Outer edge raise (meters). Default: 1, minimum: 0.0. Use period (.) for decimal separator** (Number)
- **Transversal concavity depth (meters). Default: 0.15, positive = concave. Range: 0.0 to 1.0. Use period (.) for decimal separator** (Number)
- **Back slope angle (degrees). Angle of ribbon sidewall from horizontal. Default: 45. Higher = steeper, lower = flatter. Range: 15–90. Use period (.) for decimal.** (Number)
- **Tread edge bevel size (m). Chamfers the left and right top edges. 0 = no bevel. Use period (.) for decimal.** (Number)
- **Ribbon width (meters). Default: 2.0, minimum: 0.5. Used when no field is set. Use period (.) for decimal separator** (Number)
- **Ribbon width from attribute (optional)** (Field)
- **Texture image file (optional)** (File)
- **Simple shapes mode (auto-calculate table and sigmoid lengths from heights)** (Boolean)
- **Shape distances from start (comma-separated, meters). Use period (.) for decimals. Use "end" for line end** (String)
- **Shape heights (comma-separated, meters). Use period (.) for decimals. Minimum: 0.25** (String)
- **Shape table lengths (comma-separated, meters) - flat portion between sigmoid start and end. Use period (.) for decimals** (String)
- **Sigmoid transition lengths (comma-separated, meters) - length of sigmoid curve at start and end of each shape. Must match number of shape distances. Use period (.) for decimals** (String)
- **Shape arc radius factor. In-run ramp uses a circular arc. L = in-run length; R = max(factor × L, L). 1 = quarter circle (sharpest), >1 = gentler arc. Min: 1.0. Default: 1.0.** (Number)
- **Ribbon locations GeoJSON output file (optional)** (FileDestination)
- **Export format** (Enum)

### Parameters
- **Input line layer** (VectorLayer)
- **Selected features only** (Boolean, default True)
- **DEM layer** (RasterLayer)
- **Trail raise height (meters). Default: 0.5, minimum: 0.1. Use period (.) for decimal separator** (Number, default 0.5)
- **Outer edge raise (meters). Default: 1, minimum: 0.0. Use period (.) for decimal separator** (Number, default 1)
- **Transversal concavity depth (meters). Default: 0.15, positive = concave. Range: 0.0 to 1.0. Use period (.) for decimal separator** (Number, default 0.15)
- **Back slope angle (degrees). Angle of ribbon sidewall from horizontal. Default: 45. Higher = steeper, lower = flatter. Range: 15–90. Use period (.) for decimal.** (Number, default 45.0)
- **Tread edge bevel size (m). Chamfers the left and right top edges. 0 = no bevel. Use period (.) for decimal.** (Number, default 0.0)
- **Ribbon width (meters). Default: 2.0, minimum: 0.5. Used when no field is set. Use period (.) for decimal separator** (Number, default 2)
- **Ribbon width from attribute (optional)** (Field)
- **Texture image file (optional)** (File)
- **Simple shapes mode (auto-calculate table and sigmoid lengths from heights)** (Boolean, default False)
- **Shape distances from start (comma-separated, meters). Use period (.) for decimals. Use "end" for line end** (String, default '')
- **Shape heights (comma-separated, meters). Use period (.) for decimals. Minimum: 0.25** (String, default '')
- **Shape table lengths (comma-separated, meters) - flat portion between sigmoid start and end. Use period (.) for decimals** (String, default '')
- **Sigmoid transition lengths (comma-separated, meters) - length of sigmoid curve at start and end of each shape. Must match number of shape distances. Use period (.) for decimals** (String, default '')
- **Shape arc radius factor. In-run ramp uses a circular arc. L = in-run length; R = max(factor × L, L). 1 = quarter circle (sharpest), >1 = gentler arc. Min: 1.0. Default: 1.0.** (Number, default 1.0)
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
    'TRANSVERSAL_CONCAVITY': ...,
    'BACK_SLOPE_ANGLE_DEG': ...,
    'TREAD_BEVEL_SIZE': ...,
    'RIBBON_WIDTH': ...,
    'RIBBON_WIDTH_FIELD': ...,
    'TEXTURE_FILE': ...,
    'SIMPLE_SHAPES': ...,
    'SHAPE_DISTANCES': ...,
    'SHAPE_HEIGHTS': ...,
    'SHAPE_TABLE_LENGTHS': ...,
    'SIGMOID_LENGTHS': ...,
    'SHAPE_ARC_RADIUS_FACTOR': ...,
    'OUTPUT_GEOJSON': ...,
    'EXPORT_FORMAT': ...,
    'OUTPUT': 'memory:'
})
```
