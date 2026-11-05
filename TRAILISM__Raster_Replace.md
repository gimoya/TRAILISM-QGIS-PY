## TRAILISM: Raster Replace DSM with DEM

TRAILISM ALGORITHM:\n\nReplaces DSM raster values with DEM raster values within polygon areas.\nDSM (Digital Surface Model) includes vegetation and buildings.\nDEM (Digital Elevation Model) represents the actual terrain surface.\n\nParameters:\n- DSM raster: Input DSM raster layer (base surface model with vegetation/buildings)\n- DEM raster: Input DEM raster layer (terrain values to replace DSM with)\n- Polygon layer: Defines areas where replacement occurs\n- Level ground: If checked, uses DEM value at polygon centroid for entire polygon.\n  If unchecked, uses actual DEM values for each cell within polygons.\n\nOutput: New raster with DSM values replaced by DEM values in polygon areas.

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
