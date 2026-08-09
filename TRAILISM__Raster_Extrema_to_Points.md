## TRAILISM: Raster Extrema to Points

### Algorithm id
`raster_extrema`

### Inputs
- **Input raster layer** (RasterLayer)
- **Polygon layer - restrict sampling to areas within polygons** (FeatureSource)
- **Extrema type** (Enum)

### Parameters
- **Input raster layer** (RasterLayer)
- **Polygon layer - restrict sampling to areas within polygons** (FeatureSource)
- **Extrema type** (Enum)
- **Output points** (FeatureSink)

### Outputs
- **Output points**

### Usage
- Python:
```python
import processing
processing.run('raster_extrema', {
    'INPUT_RASTER': ...,
    'POLYGONS': ...,
    'EXTREMA_TYPE': ...,
    'OUTPUT': ...,
    'OUTPUT': 'memory:'
})
```
