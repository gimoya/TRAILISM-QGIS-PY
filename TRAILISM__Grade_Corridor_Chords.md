## TRAILISM: Grade Corridor Chords

TRAILISM ALGORITHM:\n\nExact-grade downhill path from map-click start to map-click end inside a corridor.\n\nMethod:\n- Local aspect/slope from DEM sampled at cruise chord length (ignores micro relief between).\n- Grade from chord start/end DEM only: (z_a - z_b) / |ab| = target.\n- Grow on one grade diagonal (L or R) until the next step exits the corridor; clip to polygon edge.\n- Try side switch; if the switch chord is outside, backtrack one node and retry. If backtrack reaches the origin, try the other opening side.\n\nStart/end → corridor CRS; output in corridor CRS. Metres CRS required.

### Algorithm id
`grade_corridor_chords`

### Inputs
- **DEM raster** (RasterLayer)
- **Corridor polygon** (FeatureSource)
- **Start point (click on map)** (Point)
- **End point (click on map)** (Point)
- **Target grade (%)** (Number)
- **Cruise chord length (m)** (Number)
- **Minimum chord length (m)** (Number)
- **Maximum steps per attempt** (Number)
- **Output all candidate graph branches** (Boolean)

### Parameters
- **DEM raster** (RasterLayer)
- **Corridor polygon** (FeatureSource)
- **Start point (click on map)** (Point)
- **End point (click on map)** (Point)
- **Target grade (%)** (Number, default 8.0)
- **Cruise chord length (m)** (Number, default 30.0)
- **Minimum chord length (m)** (Number, default 2.0)
- **Maximum steps per attempt** (Number, default 2000)
- **Output all candidate graph branches** (Boolean, default False)
- **Output path** (FeatureSink)
- **Graph branches (one LineString per branch)** (FeatureSink)

### Outputs
- **Output path**
- **Graph branches (one LineString per branch)**

### Usage
- Python:
```python
import processing
processing.run('grade_corridor_chords', {
    'DEM': ...,
    'CORRIDOR': ...,
    'START': ...,
    'END': ...,
    'GRADE_PCT': ...,
    'STEP_LENGTH': ...,
    'MIN_CHORD_LENGTH': ...,
    'MAX_STEPS': ...,
    'KEEP_GRAPH_BRANCHES': ...,
    'OUTPUT': ...,
    'GRAPH_OUTPUT': ...,
    'OUTPUT': 'memory:'
})
```
