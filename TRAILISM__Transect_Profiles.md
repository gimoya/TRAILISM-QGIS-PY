## TRAILISM: Cross-Section Profiles

At each station along the input axis, the tool samples a DEM cross-section from s=-R (left) to s=+R (right) around the axis point (use to narrow down to a section of interest, ATTENTION: if R is too large, the tool may run out of memory, if R is too small, no intersection of sideslope with DEM is possible and C points in outputs will be missing!). The axis elevation z0 is taken from DEM plus optional dz, unless overridden by full cut/fill modes. The trail surface is flat between edges (no crossfall). From each trail edge, a side slope with angle beta_1 to horizontal is extended outward (up or down is chosen to hit the DEM first). Intersections define C points and the horizontal distance slope_run.\nParameters: STEP_M (station spacing, m); WIDTH_FIELD or WIDTH_FALLBACK (trail width, m); PROFILE_RANGE_M R (m) for sampling extent; BETA1_DEG (deg) controls side slope tan(beta_1)=rise/run; DEM_OFFSET (m) shifts DEM vertically; AXIS_DZ or AXIS_DZ_FIELD (m) adjusts z0 as DEM+dz (ignored if full cut/fill); AXIS_FULL_CUT/FILL set z0 to min/max DEM at edges; SURF_* parameters only affect SVG visualization (trapezoid body h_A, side angle to horizontal, optional top crossfall q). SVG_NAME_FIELD selects the per-feature token for Sys_Schnitt_<token>_dz±X_Y_s_<station>.svg files under SECTIONS/; effective dz per feature is included in the filename (e.g. _dz-0_75_, _dz+0_25_). 'SVG: use all stations' exports an SVG at every STEP_M station.\nOutputs: points (axis center, edges, C points) with attributes including cut/fill kind, slope_run, z, z_base, z_diff, per-interval signed areas (area_cut>0, area_fill<0). Lines: BC segments per side; chains for C (left/right) and trail edges.

### Algorithm id
`profiles_and_slope_intersections`

### Inputs
- **Input trail centerline (line)** (FeatureSource)
- **DEM raster** (RasterLayer)
- **Side slope angle beta_1 (degrees; 0–89.999). tan(beta_1) gives rise/run; legend shows nearest simple ratio.** (Number)
- **Station spacing along axis (meters; ≥ 1). Stations are multiples of this step.** (Number)
- **Width attribute field (meters). If set and valid per feature, overrides the default width.** (Field)
- **Default trail width (meters; ≥ 0). Used when width field is missing/invalid.** (Number)
- **Profile half-range R (meters; ≥ 1). DEM is sampled from -R (left) to +R (right) around axis.** (Number)
- **DEM vertical offset for top-soil removal (meters; -1 to 0). Applied globally to DEM for all calculations and plots.** (Number)
- **Default axis elevation offset dz (meters; -2.5 to 2.5). z0 := DEM + dz unless overridden by full cut/fill mode or dz field.** (Number)
- **Axis elevation offset field dz (meters). If set and valid per feature, overrides default dz.** (Field)
- **Axis mode: full cut — set z0 to the lower DEM edge (min). Ignores dz.** (Boolean)
- **Axis mode: full fill — set z0 to the higher DEM edge (max). Ignores dz.** (Boolean)
- **Surface thickness height h_A (meters; 0–1). Visual overlay only; does not affect cut/fill.** (Number)
- **Surface thickness field (meters). If present per feature, overrides h_A.** (Field)
- **Side face angle to horizontal (degrees; 20–90). 90° = vertical sides (rectangle). Visual overlay only.** (Number)
- **Surface top crossfall (visual only)** (Enum)
- **Surface top crossfall tilt q (degrees; 0–7). Used only when crossfall ≠ Flat.** (Number)
- **SVG filename field — unique value per line feature. Token in Sys_Schnitt_<value>_dz+0_25_s_<station>.svg (dz from AXIS_DZ / AXIS_DZ_FIELD). Pick a field with distinct values per segment; duplicate values overwrite earlier SVGs (warning in log). Falls back to feat_<id> if unset, missing, or empty.** (Field)
- **SVG: use all stations** (Boolean)
- **SVG stations to plot (meters, comma-separated; exact match). Capped to each line part length [0, L].** (String)
- **SVG scale: pixels per meter (30–300). Controls both axes; aspect is 1:1.** (Number)

### Parameters
- **Input trail centerline (line)** (FeatureSource)
- **DEM raster** (RasterLayer)
- **Side slope angle beta_1 (degrees; 0–89.999). tan(beta_1) gives rise/run; legend shows nearest simple ratio.** (Number, default 33.7)
- **Station spacing along axis (meters; ≥ 1). Stations are multiples of this step.** (Number, default 10.0)
- **Width attribute field (meters). If set and valid per feature, overrides the default width.** (Field)
- **Default trail width (meters; ≥ 0). Used when width field is missing/invalid.** (Number, default 2.0)
- **Profile half-range R (meters; ≥ 1). DEM is sampled from -R (left) to +R (right) around axis.** (Number, default 2.5)
- **DEM vertical offset for top-soil removal (meters; -1 to 0). Applied globally to DEM for all calculations and plots.** (Number, default -0.15)
- **Default axis elevation offset dz (meters; -2.5 to 2.5). z0 := DEM + dz unless overridden by full cut/fill mode or dz field.** (Number, default -0.1)
- **Axis elevation offset field dz (meters). If set and valid per feature, overrides default dz.** (Field)
- **Axis mode: full cut — set z0 to the lower DEM edge (min). Ignores dz.** (Boolean, default False)
- **Axis mode: full fill — set z0 to the higher DEM edge (max). Ignores dz.** (Boolean, default False)
- **Surface thickness height h_A (meters; 0–1). Visual overlay only; does not affect cut/fill.** (Number, default 0.2)
- **Surface thickness field (meters). If present per feature, overrides h_A.** (Field)
- **Side face angle to horizontal (degrees; 20–90). 90° = vertical sides (rectangle). Visual overlay only.** (Number, default 45.0)
- **Surface top crossfall (visual only)** (Enum)
- **Surface top crossfall tilt q (degrees; 0–7). Used only when crossfall ≠ Flat.** (Number, default 3)
- **SVG filename field — unique value per line feature. Token in Sys_Schnitt_<value>_dz+0_25_s_<station>.svg (dz from AXIS_DZ / AXIS_DZ_FIELD). Pick a field with distinct values per segment; duplicate values overwrite earlier SVGs (warning in log). Falls back to feat_<id> if unset, missing, or empty.** (Field)
- **SVG: use all stations** (Boolean, default False)
- **SVG stations to plot (meters, comma-separated; exact match). Capped to each line part length [0, L].** (String, default "")
- **SVG scale: pixels per meter (30–300). Controls both axes; aspect is 1:1.** (Number, default 180.0)
- **Points: axis center, trail edges, slope toes/heels (with areas)** (FeatureSink, default QgsProcessing.TEMPORARY_OUTPUT)
- **BC segments per side (line)** (FeatureSink, default QgsProcessing.TEMPORARY_OUTPUT)
- **Chains: trail edges and C (line)** (FeatureSink, default QgsProcessing.TEMPORARY_OUTPUT)

### Outputs
- **Points: axis center, trail edges, slope toes/heels (with areas)**
- **BC segments per side (line)**
- **Chains: trail edges and C (line)**

### Usage
- Python:
```python
import processing
processing.run('profiles_and_slope_intersections', {
    'LINE_LAYER': ...,
    'DEM': ...,
    'BETA1_DEG': ...,
    'STEP_M': ...,
    'WIDTH_FIELD': ...,
    'WIDTH_FALLBACK': ...,
    'PROFILE_RANGE_M': ...,
    'DEM_OFFSET': ...,
    'AXIS_DZ': ...,
    'AXIS_DZ_FIELD': ...,
    'AXIS_FULL_CUT': ...,
    'AXIS_FULL_FILL': ...,
    'SURF_HEIGHT': ...,
    'SURF_HEIGHT_FIELD': ...,
    'SURF_CROSSFALL_DEG': ...,
    'SURF_TOP_CROSSFALL_MODE': ...,
    'SURF_TOP_TILT_DEG': ...,
    'SVG_NAME_FIELD': ...,
    'SVG_ALL_STATIONS': ...,
    'SVG_STATIONS': ...,
    'SVG_PX_PER_M': ...,
    'OUT_POINTS': ...,
    'OUT_BC_SEGMENTS': ...,
    'OUT_CHAIN_C': ...,
    'OUTPUT': 'memory:'
})
```
