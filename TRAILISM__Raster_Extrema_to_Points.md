## TRAILISM: Raster Extrema to Points

TRAILISM ALGORITHM:\n\nFinds minimum and maximum raster values within each polygon feature.\nOutputs points at the locations of these extrema.\n\nParameters:\n- Extrema type: find maxima (peaks), minima (depressions), or both\n\nPolygon layer is required - finds extrema within each polygon feature.

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
