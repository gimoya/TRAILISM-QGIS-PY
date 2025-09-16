## TRAILISM: Cross-Section Profiles

Samples DEM cross-sections perpendicular to a centerline at fixed stations, finds slope intersections from each trail edge (C points), computes per-station cut/fill areas, and optionally exports equal-scale SVGs.

### Algorithm id
`profiles_and_slope_intersections` (check in your Processing pane)

### Inputs
- **LINE_LAYER**: Trail centerline (line)
- **DEM**: DEM raster layer

### Core parameters
- **STEP_M** (Double, default 10): Station spacing along axis.
- **WIDTH_FIELD** (Numeric field, optional): Per-feature width override.
- **WIDTH_FALLBACK** (Double, default 2.0): Default width when field missing.
- **PROFILE_RANGE_M** (Double, default 2.5): Half-range R for sampling (left/right extent).
- **BETA1_DEG** (Double, default 33.7): Side slope angle to horizontal; tan(β1) = rise/run.
- **DEM_OFFSET** (Double, default -0.15): Global vertical shift applied to DEM for all calcs and plots.
- **AXIS_DZ** (Double, default -0.1): Additional z-offset at axis center (ignored in full cut/fill modes).
- **AXIS_FULL_CUT / AXIS_FULL_FILL** (Bool): Override z0 to min/max DEM at edges, ignoring dz.

### SVG (optional, visualization only)
- **SURF_HEIGHT** (Double, default 0.2): Surface thickness for the trapezoid body (visual only).
- **SURF_HEIGHT_FIELD** (Numeric field, optional): Per-feature override of thickness.
- **SURF_CROSSFALL_DEG** (Double, default 45): Side face angle to horizontal; 90° = vertical.
- **SURF_TOP_CROSSFALL_MODE** (Enum: Flat | Inward | Outward; default Outward)
- **SURF_TOP_TILT_DEG** (Int 0–7, default 3): Tilt magnitude for the top when crossfall ≠ Flat.
- **SVG_OPEN** (Bool, default True): Open exported SVGs.
- **SVG_ALL_STATIONS** (Bool, default False): Export all stations (capped at 100 SVGs/run).
- **SVG_STATIONS** (String): Comma-separated list of exact station meters to plot.
- **SVG_PX_PER_M** (Double, default 180): 1:1 scaling in pixels per meter.

### Outputs
- **OUT_POINTS** (Point): Axis center, edges, and C points. Fields include station, type, side, cut/fill, slope_run, z, z_base, z_diff, svg_sel, and total per-station areas (area_cut_total, area_fill_total).
- **OUT_BC_SEGMENTS** (Line): Short BC segments from edges to C; also transect lines for plotted stations.
- **OUT_CHAIN_C** (Line): Chains of C points (left/right) and trail edge chains; `OUT_CHAIN_EDGES` is the same layer for compatibility.
- **OUT_AREAS**: Alias of points layer (kept for compatibility).

### Usage
- GUI: Processing Toolbox → “TRAILISM: Cross-Section Profiles”.
- Python:
```python
import processing
processing.run('profiles_and_slope_intersections', {
    'LINE_LAYER': axis_layer,
    'DEM': dem_layer,
    'STEP_M': 10.0,
    'WIDTH_FIELD': None,           # or 'width'
    'WIDTH_FALLBACK': 2.0,
    'PROFILE_RANGE_M': 2.5,
    'BETA1_DEG': 33.7,
    'DEM_OFFSET': -0.15,
    'AXIS_DZ': -0.1,
    'AXIS_FULL_CUT': False,
    'AXIS_FULL_FILL': False,
    'SURF_HEIGHT': 0.2,
    'SURF_HEIGHT_FIELD': None,
    'SURF_CROSSFALL_DEG': 45.0,
    'SURF_TOP_CROSSFALL_MODE': 2,  # 0 Flat, 1 Inward, 2 Outward
    'SURF_TOP_TILT_DEG': 3,
    'SVG_OPEN': True,
    'SVG_ALL_STATIONS': False,
    'SVG_STATIONS': '',
    'SVG_PX_PER_M': 180.0,
    'OUT_POINTS': 'memory:',
    'OUT_BC_SEGMENTS': 'memory:',
    'OUT_CHAIN_C': 'memory:'
})
```

### Tips
- Keep R modest; very large R increases memory and can miss intersections at realistic side slopes.
- If C is missing on one side, that side’s chain will fall back to the edge at that station to keep polylines continuous.
- DEM/axis CRS mismatch is handled internally via coordinate transform.

### Styling
- Outputs are pre-styled white and small for overlay on maps; points include a rule to highlight SVG-selected stations.


