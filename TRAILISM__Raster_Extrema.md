## TRAILISM: Raster Extrema Points

Finds minimum and maximum raster values within each polygon feature. Outputs points at the locations of these extrema.

### Algorithm id
`raster_extrema` (check your Processing id via tooltip)

### Inputs
- **INPUT_RASTER**: Raster layer
- **POLYGONS**: Polygon layer - restrict sampling to areas within polygons (required)

### Parameters
- **EXTREMA_TYPE** (Enum, default 0): Type of extrema to find
  - 0: Maxima only (peaks)
  - 1: Minima only (depressions)
  - 2: Both

### Outputs
- **OUTPUT**: Point layer with extrema locations. Fields include:
  - `extrema_type`: "maximum" or "minimum"
  - `value`: The raster value at the extrema point
  - `polygon_id`: ID of the polygon feature containing the extrema

### Behavior
- Processes each polygon feature separately to find extrema within its bounds.
- Finds all raster cells within each polygon and identifies minimum/maximum values.
- When multiple cells share the same extremum value, outputs a point at the centroid of those cells.
- Automatically handles coordinate reference system transformations between raster and polygon layers.
- Output points are styled with rule-based rendering: red circles for maxima, green circles for minima.

### Usage
- GUI: Processing Toolbox → "TRAILISM: Raster Extrema Points".
- Python:
```python
import processing
processing.run('raster_extrema', {
    'INPUT_RASTER': raster_layer,
    'POLYGONS': polygon_layer,
    'EXTREMA_TYPE': 2,  # 0 = Maxima only, 1 = Minima only, 2 = Both
    'OUTPUT': 'memory:'
})
```

### Tips
- Polygon layer is required - the algorithm finds extrema within each polygon feature separately.
- Use this to identify peaks and depressions within specific areas of interest (e.g., watershed boundaries, administrative units).
- If multiple cells share the same extremum value, the output point represents the centroid of those cells.
- The algorithm automatically handles CRS mismatches between raster and polygon layers.

