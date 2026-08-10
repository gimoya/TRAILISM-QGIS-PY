## TRAILISM: Grade Corridor Chords

TRAILISM ALGORITHM:\n\nWalk downhill inside a corridor:\n1. Full circle at chord length → steepest-down orientation (ties → angular mid).\n2. UI L/R semicircle (0…180°) → closest |grade − target|.\n3. Ray inside → next node; ray exits → edge intersection, flip side, continue.\n4. Two consecutive short chords / no candidate → drop tip, rewind, opposite exit side (once per node/side).\n5. End when within UI buffer of end point, or rewind exhausts back at start.\n\nOutput: node points + all try paths as LineStrings + debug fall-line ray and 180° pick fan (best path).

### Algorithm id
`grade_corridor_chords`

### Inputs
- **DEM raster** (RasterLayer)
- **Corridor polygon** (FeatureSource)
- **Start point (click on map)** (Point)
- **End point (click on map)** (Point)
- **End reach buffer (m)** (Number)
- **Chord length (m)** (Number)
- **Side of steepest-down orientation** (Enum)
- **Desired grade (%)** (Number)

### Parameters
- **DEM raster** (RasterLayer)
- **Corridor polygon** (FeatureSource)
- **Start point (click on map)** (Point)
- **End point (click on map)** (Point)
- **End reach buffer (m)** (Number, default 5.0)
- **Chord length (m)** (Number, default 30.0)
- **Side of steepest-down orientation** (Enum)
- **Desired grade (%)** (Number, default 8.0)
- **Nodes** (FeatureSink)
- **Path tries (all LineStrings)** (FeatureSink)
- **Debug: fall-line + 180° pick fan** (FeatureSink)

### Outputs
- **Nodes**
- **Path tries (all LineStrings)**
- **Debug: fall-line + 180° pick fan**

### Usage
- Python:
```python
import processing
processing.run('grade_corridor_chords', {
    'DEM': ...,
    'CORRIDOR': ...,
    'START': ...,
    'END': ...,
    'END_BUFFER': ...,
    'STEP_LENGTH': ...,
    'SIDE': ...,
    'GRADE_PCT': ...,
    'OUTPUT': ...,
    'LINE_OUTPUT': ...,
    'FALLLINE_OUTPUT': ...,
    'OUTPUT': 'memory:'
})
```
