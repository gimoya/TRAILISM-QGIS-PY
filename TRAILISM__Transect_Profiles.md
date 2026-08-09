## TRAILISM__Transect_Profiles

### Inputs
- **Input trail centerline (line)** (FeatureSource)
- **DEM raster** (RasterLayer)
- **Side slope angle beta_1 (degrees; 0–89.999). tan(beta_1) gives rise/run; legend shows nearest simple ratio.** (Number)
- **Station spacing along axis (meters; ≥ 1). Stations are multiples of this step.** (Number)
- **Width attribute field (meters). If set and valid per feature, overrides the default width.** (Field)
- **Default trail width (meters; ≥ 0). Used when width field is missing/invalid.** (Number)
- **Profile half-range R (meters; ≥ 1). DEM is sampled from -R (left) to +R (right) around axis.** (Number)
- **DEM vertical offset (meters; -1 to 0). Applied globally to DEM for all calculations and plots.** (Number)
- **Default axis elevation offset dz (meters; -2.5 to 2.5). z0 := DEM + dz unless overridden by full cut/fill mode or dz field.** (Number)
- **Axis elevation offset field dz (meters). If set and valid per feature, overrides default dz.** (Field)
- **Axis mode: full cut — set z0 to the lower DEM edge (min). Ignores dz.** (Boolean)
- **Axis mode: full fill — set z0 to the higher DEM edge (max). Ignores dz.** (Boolean)
- **Surface thickness height h_A (meters; 0–1). Visual overlay only; does not affect cut/fill.** (Number)
- **Surface thickness field (meters). If present per feature, overrides h_A.** (Field)
- **Side face angle to horizontal (degrees; 20–90). 90° = vertical sides (rectangle). Visual overlay only.** (Number)
- **Surface top crossfall (visual only)** (Enum)
- **Surface top crossfall tilt q (degrees; 0–7). Used only when crossfall ≠ Flat.** (Number)
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
- **DEM vertical offset (meters; -1 to 0). Applied globally to DEM for all calculations and plots.** (Number, default -0.15)
- **Default axis elevation offset dz (meters; -2.5 to 2.5). z0 := DEM + dz unless overridden by full cut/fill mode or dz field.** (Number, default -0.1)
- **Axis elevation offset field dz (meters). If set and valid per feature, overrides default dz.** (Field)
- **Axis mode: full cut — set z0 to the lower DEM edge (min). Ignores dz.** (Boolean, default False)
- **Axis mode: full fill — set z0 to the higher DEM edge (max). Ignores dz.** (Boolean, default False)
- **Surface thickness height h_A (meters; 0–1). Visual overlay only; does not affect cut/fill.** (Number, default 0.2)
- **Surface thickness field (meters). If present per feature, overrides h_A.** (Field)
- **Side face angle to horizontal (degrees; 20–90). 90° = vertical sides (rectangle). Visual overlay only.** (Number, default 45.0)
- **Surface top crossfall (visual only)** (Enum)
- **Surface top crossfall tilt q (degrees; 0–7). Used only when crossfall ≠ Flat.** (Number, default 3)
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
- **SVG cross-section diagrams** (optional): written to `SECTIONS/` next to the QGIS project file
