## TRAILISM: Raster Replace DSM with DEM

Replaces DSM raster values with DEM raster values within polygon areas. DSM (Digital Surface Model) includes vegetation and buildings. DEM (Digital Elevation Model) represents the actual terrain surface. Optionally levels the ground by using a single DEM value (sampled at polygon centroid) for all cells within each polygon.

### Algorithm id
`raster_replace` (check your Processing id via tooltip)

### Inputs
- **DSM_RASTER**: DSM raster layer (Digital Surface Model - base raster with vegetation/buildings to be modified)
- **DEM_RASTER**: DEM raster layer (Digital Elevation Model - source of terrain replacement values)
- **POLYGONS**: Polygon layer - defines areas where replacement occurs (required)

### Parameters
- **LEVEL_GROUND** (Boolean, default False): Level ground option
  - If **unchecked** (False): Uses actual DEM values for each cell within polygons
  - If **checked** (True): Uses DEM value at polygon centroid for entire polygon area (levels the ground)

### Outputs
- **OUTPUT**: Raster layer with DSM values replaced by DEM values in polygon areas

### Behavior
- Processes each polygon feature separately
- Copies DSM raster to output first
- For each polygon:
  - If **level_ground** is False: Replaces DSM values with corresponding DEM values for each cell within the polygon
  - If **level_ground** is True: Samples DEM at polygon centroid and uses that single value for all cells within the polygon
- Automatically handles coordinate reference system transformations between rasters and polygon layer
- Preserves DSM data type and nodata values in output
- Output raster uses LZW compression and tiled format for efficient storage

### Usage
- GUI: Processing Toolbox → "TRAILISM: Raster Replace DSM with DEM".
- Python:
```python
import processing
processing.run('raster_replace', {
    'DSM_RASTER': dsm_raster_layer,
    'DEM_RASTER': dem_raster_layer,
    'POLYGONS': polygon_layer,
    'LEVEL_GROUND': False,  # False = use actual DEM values, True = level ground
    'OUTPUT': 'path/to/output.tif'
})
```

### Tips
- Polygon layer is required - defines areas where replacement occurs
- Use **level_ground = False** when you want to preserve DEM terrain detail within polygons
- Use **level_ground = True** when you want to create flat/level areas (e.g., for construction sites, platforms)
- Both DSM and DEM rasters should have the same cell size and alignment for best results
- The algorithm automatically handles CRS mismatches between rasters and polygon layer
- Output preserves the data type of the input DSM raster

