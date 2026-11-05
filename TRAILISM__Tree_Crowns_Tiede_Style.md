## TRAILISM: Tree peaks (Tiede local maxima)

Tree top localization using Tiede-style local maxima: CD = min(a + b*TH, cap), search radius R = CD/2 at each focal cell (coarse max filter, then per-cell disk; ties: row-major).\n\nEdge padding (meters) extends the grid with low synthetic cells so disks near the raster border are fully sampled; outputs are only written for the original extent.\n\nNo crown delineation — peak points only. Defaults a=1.54, b=0.123 from Tiede et al. (2005).\nInput: DSM or CHM/nDOM; optional DTM for nDOM = DSM − DTM; optional extent.\nRequires scipy.ndimage for the coarse max filter.

### Algorithm id
`tree_peaks_tiede_style`

### Inputs
- **Input raster (DSM or canopy height CHM/nDOM)** (RasterLayer)
- **DTM (optional) — if set, nDOM = DSM − DTM** (RasterLayer)
- **Extent (optional) — limit processing to this area** (Extent)
- **Crown model: a in CD = a + b·TH (meters)** (Number)
- **Crown model: b in CD = a + b·TH** (Number)
- **Maximum crown diameter cap (meters)** (Number)
- **Minimum canopy height — preliminary mask (meters)** (Number)
- **Edge padding (meters) — pad grid for full disk sampling** (Number)
- **Output tree peaks (local maxima, row-major tie-break)** (VectorDestination)

### Parameters
- **Input raster (DSM or canopy height CHM/nDOM)** (RasterLayer)
- **DTM (optional) — if set, nDOM = DSM − DTM** (RasterLayer)
- **Extent (optional) — limit processing to this area** (Extent)
- **Crown model: a in CD = a + b·TH (meters)** (Number, default 1.54)
- **Crown model: b in CD = a + b·TH** (Number, default 0.123)
- **Maximum crown diameter cap (meters)** (Number, default 10.0)
- **Minimum canopy height — preliminary mask (meters)** (Number, default 5.0)
- **Edge padding (meters) — pad grid for full disk sampling** (Number, default 10.0)
- **Output tree peaks (local maxima, row-major tie-break)** (VectorDestination)

### Usage
- Python:
```python
import processing
processing.run('tree_peaks_tiede_style', {
    'INPUT_RASTER': ...,
    'INPUT_DTM': ...,
    'EXTENT': ...,
    'COEFF_A': ...,
    'COEFF_B': ...,
    'MAX_CROWN_DIAMETER_M': ...,
    'MIN_TREE_HEIGHT_M': ...,
    'EDGE_PAD_M': ...,
    'OUTPUT_PEAKS': ...,
    'OUTPUT': 'memory:'
})
```
