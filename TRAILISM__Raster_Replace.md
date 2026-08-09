## TRAILISM: Raster Replace DSM with DEM

### Algorithm id
`raster_replace`

### Inputs
- **DSM raster layer (Digital Surface Model)** (RasterLayer)
- **DEM raster layer (Digital Elevation Model)** (RasterLayer)
- **Polygon layer - defines areas for replacement** (FeatureSource)
- **Level ground - use DEM value at polygon centroid for entire polygon** (Boolean)
- **Output raster** (RasterDestination)

### Parameters
- **DSM raster layer (Digital Surface Model)** (RasterLayer)
- **DEM raster layer (Digital Elevation Model)** (RasterLayer)
- **Polygon layer - defines areas for replacement** (FeatureSource)
- **Level ground - use DEM value at polygon centroid for entire polygon** (Boolean, default False)
- **Output raster** (RasterDestination)

### Usage
- Python:
```python
import processing
processing.run('raster_replace', {
    'DSM_RASTER': ...,
    'DEM_RASTER': ...,
    'POLYGONS': ...,
    'LEVEL_GROUND': ...,
    'OUTPUT': ...,
    'OUTPUT': 'memory:'
})
```
