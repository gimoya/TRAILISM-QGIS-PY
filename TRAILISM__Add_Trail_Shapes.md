## TRAILISM: Add Trail Shapes

Creates 3D trail models with curve banking, transversal concavity, and optional sigmoid "shape" bumps. Exports ribbon meshes in GLTF/GLB format for 3D visualization in Three.js or other 3D viewers. All Z values are derived from DEM sampling, completely ignoring any Z component from the input linestring geometry.

### Algorithm id
`add_trail_shapes` (check your Processing id via tooltip)

### Inputs
- **INPUT**: Line layer (trail centerline)
- **DEM**: DEM raster layer (required) - used for sampling base Z values

### Parameters

#### Basic Parameters
- **ONLY_RAISE_TRAIL** (Boolean, default False): If checked, skip shapes and curve banking, just raise trail above terrain
- **SELECTED_FEATURES_ONLY** (Boolean, default True): Process only selected features
- **Z_RAISE** (Double, default 0.5): Trail raise height in meters above terrain. Minimum: 0.0
- **RIBBON_WIDTH** (Double, default 3.0): Base width of the exported ribbon mesh in meters. Minimum: 0.5
- **DENSIFY_DISTANCE** (Double, default 0.25): Maximum distance between vertices in meters. Minimum: 0.05
- **EXPORT_FORMAT** (Enum, default GLB): Export format - GLTF or GLB
- **TEXTURE_FILE** (File, optional): Texture image file (PNG, JPG, JPEG) for the ribbon mesh

#### Curve Banking Parameters
- **OUTER_EDGE_RAISE** (Double, default 1.0, optional): Height to raise the outer edge of curves in meters. Minimum: 0.0
- **TRANSVERSAL_CONCAVITY** (Double, default 0.25, optional): Transversal concavity depth in meters. Positive = concave. Range: 0.0 to 1.0

**Note**: Curve detection uses internal constants (window distance: 15.0m, threshold angle: 10.0°) and is not exposed as parameters.

#### Shape Parameters
- **SIMPLE_SHAPES** (Boolean, default False): Simple shapes mode
  - If **checked**: Only distances and heights are required. Table and sigmoid lengths are auto-calculated:
    - Table length = height × 3 - 1.5
    - Sigmoid length = height × 4
    - Minimum height: 0.5 meters
  - If **unchecked**: All shape parameters must be provided manually

- **SHAPE_DISTANCES** (String, optional): Comma-separated distances from start in meters
  - Use period (.) for decimal values
  - Use "end" for line end, "0" for start
  - Example: "10,30,end"

- **SHAPE_HEIGHTS** (String, optional): Comma-separated heights for each shape in meters
  - Must match number of distances
  - Minimum: 0.5 meters (enforced in simple mode)
  - Use period (.) for decimal values
  - Example: "1.5,2.0,1.0"

- **SHAPE_TABLE_LENGTHS** (String, optional): Comma-separated flat portion lengths in meters
  - Required only if Simple shapes mode is unchecked
  - Must match number of distances
  - Example: "3.0,4.5,3.0"

- **SIGMOID_LENGTHS** (String, optional): Comma-separated lengths for sigmoid transition curves in meters
  - Required only if Simple shapes mode is unchecked
  - Must match number of distances
  - Example: "6.0,8.0,4.0"

- **SIGMOID_STEEPNESS_K** (Double, default 4.5, optional): Controls transition steepness
  - Higher = steeper transition, lower = flatter transition
  - Range: 1.0 to 30.0

### Outputs
- **GLTF/GLB files**: Exported ribbon meshes saved to the project home directory
  - Filename format: `{layer_name}_{feature_id}.gltf` or `.glb`
  - Original layer is not modified

### Behavior

#### Z Value Handling
- **All Z values are derived from DEM sampling** - any Z component in the input linestring geometry is completely ignored
- Only XY coordinates from the input line are used
- Base Z = DEM sampled value + Z_RAISE + shape offsets

#### Adaptive Densification
- **Curves**: Vertices at 0.25m spacing (or densify_distance if smaller)
- **Straights**: Vertices at 1.25m spacing (5× densify_distance)
- Automatically detects curves and adjusts vertex density accordingly

#### Curve Banking
- Detects curves using bisector angle analysis (internal: 15m window, 10° threshold)
- Raises outer edge of curves with smooth parabolic arc
- Creates natural banking effect for trail turns

#### Transversal Concavity
- Creates concave profile across trail width
- Positive values create a "dished" trail surface
- Applied uniformly along the trail

#### Shapes (Sigmoid Bumps)
- Smooth sigmoid-shaped bumps at specified distances
- Each shape has:
  - Flat table portion (middle)
  - Sigmoid transitions at start and end
- Can be combined with curve banking for complex profiles

#### Vertex Duplication
- Duplicates vertices at sharp edges to preserve normals
- Normal threshold: 0.8 (dot product)
- Prevents unwanted smoothing at edges

#### Coordinate Precision
- Uses original geometry directly for XY positions
- Avoids precision loss from intermediate geometry reconstruction
- References original line geometry for all coordinate calculations

### Usage
- **GUI**: Processing Toolbox → "TRAILISM: Add Trail Shapes"
- **Python**:
```python
import processing
processing.run('add_trail_shapes', {
    'INPUT': line_layer,
    'DEM': dem_raster_layer,
    'ONLY_RAISE_TRAIL': False,
    'SELECTED_FEATURES_ONLY': True,
    'Z_RAISE': 0.5,
    'OUTER_EDGE_RAISE': 1.0,
    'TRANSVERSAL_CONCAVITY': 0.25,
    'RIBBON_WIDTH': 3.0,
    'SIMPLE_SHAPES': True,
    'SHAPE_DISTANCES': '10,30,end',
    'SHAPE_HEIGHTS': '1.5,2.0,1.0',
    'DENSIFY_DISTANCE': 0.25,
    'EXPORT_FORMAT': 1,  # 0 = GLTF, 1 = GLB
    'TEXTURE_FILE': None,
    'OUTPUT': None  # Files saved to project home
})
```

### Tips

#### Decimal Separator
- **IMPORTANT**: Always use period (.) as decimal separator
- Comma (,) is used only to separate list items
- GUI may display defaults with comma due to system locale, but always enter values using period

#### DEM Requirements
- DEM layer is required - all elevation values come from DEM sampling
- DEM should cover the entire trail extent
- Input linestring Z values are completely ignored

#### Simple Shapes Mode
- Recommended for most use cases
- Automatically calculates table and sigmoid lengths from heights
- Formula: table = h × 3 - 1.5, sigmoid = h × 4
- Minimum height: 0.5 meters

#### Curve Banking
- Curve detection is automatic (internal parameters)
- Outer edge raise creates natural banking on turns
- Works best with densified geometry (default 0.25m)

#### Export Format
- **GLB**: Binary format, smaller file size, recommended
- **GLTF**: Text format, human-readable JSON
- Files are saved to project home directory
- Compatible with Three.js and most 3D viewers

#### Performance
- Adaptive densification reduces vertex count in straights
- Processing time depends on trail length and densification distance
- Large trails with many shapes may take longer to process

#### Three.js Integration
- Exported models use absolute coordinates
- Compatible with Qgis2threejs plugin scenes
- Models can be loaded directly into Three.js scenes

### Example Use Cases
- Create 3D trail models for visualization
- Generate trail meshes with realistic banking on curves
- Add smooth bumps (jumps, rollers) to trail profiles
- Export trail geometry for web-based 3D viewers
- Combine multiple effects (banking + shapes + concavity) for complex trails

### Technical Notes
- Uses adaptive densification: 0.25m in curves, 1.25m in straights
- Vertex duplication threshold: 0.8 (normal dot product)
- Curve detection: 15m window, 10° threshold (internal constants)
- All Z values from DEM sampling (input geometry Z ignored)
- Original geometry referenced directly for XY precision
- GLTF/GLB export uses 32-bit indices for large meshes
- Texture coordinates automatically generated if texture file provided

