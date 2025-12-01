## TRAILISM: Raster Find Tree Tips

Creates an average slope angle raster map for tree canopy detection. For each cell, calculates the average slope angle in 8 directions by stepping through all cells up to the looking distance. Outputs the average slope angle in degrees if above minimum threshold. Positive values indicate peaks (downhill in all directions), negative values indicate depressions or slopes.

### Algorithm id
`raster_find_tree_tips` (check your Processing id via tooltip)

### Inputs
- **INPUT_RASTER**: Raster layer (typically a DEM/elevation raster)

### Parameters
- **EXTENT** (Extent, optional): Select on canvas or use current canvas view to restrict processing area. Automatically transformed to raster CRS.
- **LOOKING_DISTANCE** (Double, default 12.0): Distance to check in each direction (meters). Should be 1.5-2x typical crown diameter to ensure comparison to ground. For coniferous trees (5-12m crowns): use 10-20m. Controls how far to look when calculating average slope.
- **MIN_THRESHOLD** (Double, default 15.0): Minimum average slope angle in degrees. Only cells with average slope angle strictly above this threshold are output. Typical values: 10-15° (small trees), 15-25° (medium), 25-35° (large). Default: 15° (filters small peaks, detects medium+ trees). Range: 0.0-90.0°.

### Outputs
- **OUTPUT**: Raster layer with average slope angles in degrees. Each cell contains:
  - **Positive values (degrees)**: Peaks (downhill in all 8 directions). Higher values = steeper/more prominent peaks. Only output if above minimum threshold.
  - **NoData**: Cells below or equal to threshold, border cells (within looking distance of raster edge), or cells with incomplete paths.
  
  **Automatic Styling**: The output raster is automatically styled with:
  - Positive values (peaks) displayed in solid red
  - Zero and negative values displayed as transparent

### Behavior
- For each cell, calculates average slope angle in 8 directions (N, NE, E, SE, S, SW, W, NW)
- Steps through all cells up to the looking distance in each direction
- Collects elevation values along each path (excluding nodata cells from average)
- Calculates slope as: `(cell_elevation - average_path_elevation) / euclidean_distance`
- Converts slope to angle in degrees: `degrees(atan(slope))`
- Outputs the average of the 8 slope angles in degrees, but only if strictly above minimum threshold
- Uses Euclidean distance for diagonal directions (accurate distance calculation)
- Cells within looking distance of raster border are excluded from processing (set to NoData)
- Border cells are still used for neighbor lookups when calculating slopes for non-border cells

### Algorithm Details
- **Single-pass algorithm**: Processes each cell once
- **Efficient border handling**: Border cells excluded during processing but available for neighbor calculations
- **Pre-calculated distances**: Euclidean distances computed once per direction (optimization)
- **Coordinate conversion**: Only performed when extent filtering is used
- **Progress reporting**: Adaptive based on number of cells to process
- **Threshold filtering**: Only outputs cells with average slope angle > MIN_THRESHOLD

### Interpretation of Output Values
- **High positive values (e.g., 25-45°)**: Very prominent peaks with steep slopes in all directions
- **Low positive values (e.g., 15-25°)**: Gentle peaks or subtle local maxima
- **NoData**: Cells below threshold, border cells, or cells where paths hit raster boundaries

### Usage
- **GUI**: Processing Toolbox → "TRAILISM: Raster Find Tree Tips"
- **Python**:
```python
import processing
processing.run('raster_find_tree_tips', {
    'INPUT_RASTER': raster_layer,
    'EXTENT': None,  # Optional: extent rectangle or None
    'LOOKING_DISTANCE': 12.0,  # meters (default: 12.0)
    'MIN_THRESHOLD': 15.0,  # degrees (default: 15.0)
    'OUTPUT': '/path/to/output.tif'
})
```

### Tips
- **Looking distance**: Controls how far to check in each direction. Should be 1.5-2x typical crown diameter. Larger values find more prominent features but may miss smaller peaks. Also affects minimum spacing between detectable peaks. For coniferous trees (5-12m crowns): use 10-20m.
- **Minimum threshold**: Filters out small peaks and noise. Lower values detect more features but may include false positives. Higher values detect only prominent peaks. Typical: 10-15° (small trees), 15-25° (medium), 25-35° (large).
- **Extent parameter**: Use to restrict processing to a specific area, useful for large rasters or focused analysis.
- **Border cells**: Cells within looking distance of raster edge are automatically excluded. Ensure your raster has sufficient buffer around areas of interest.
- **Cell resolution**: For 0.5m resolution rasters with 12m looking distance, peaks can be as close as the cell resolution (0.5-1m) if they're at different elevations.
- **Visualization**: The output raster is automatically styled with positive values (peaks) in red and zero/negative values transparent. You can modify the styling in QGIS layer properties if needed.

### Example Use Cases
- Identify tree tops in high-resolution forest DEMs
- Find local peaks in terrain analysis
- Detect prominent features in elevation models
- Create slope prominence maps for trail planning
- Analyze terrain convexity/concavity

### Technical Notes
- Output values are in degrees (slope angle)
- Uses Float32 data type for output raster
- Statistics are computed automatically for proper QGIS display
- Algorithm handles nodata cells by excluding them from path averages but continuing the path
- Euclidean distance ensures accurate slope calculations for diagonal directions
- Output raster is automatically styled with a color ramp: positive values (peaks) = red, zero/negative = transparent
- Only cells with average slope angle strictly above MIN_THRESHOLD are output (others set to NoData)

