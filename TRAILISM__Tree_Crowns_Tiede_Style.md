## TRAILISM: Tree peaks (Tiede-style local maxima)

Raster analogue of the **tree-top localization** step from **Tiede, Hochleitner, Blaschke (2005)** — *A fully GIS-based workflow for tree identification and tree crown delineation using laser scanning* ([Research Studio PDF](https://www.researchstudio.at/wp-content/uploads/2020/05/Tiede-D.-Hochleitner-G.-Blaschke-T.-2005-A-fully-GIS-Based-workflow-for-tree-identification-and-tree-crown-delineation-using-laser-scanning..pdf)). Defaults **a = 1.54**, **b = 0.123** in **CD = a + b·TH** (crown diameter vs tree height); per cell **r_search_m = CD/2** after capping CD. **No crown delineation** — strict tree-peak points only (coarse filter is internal).

**Dependency:** **SciPy** (`scipy.ndimage`) is required for the coarse focal-maximum filter. Other TRAILISM scripts do not use SciPy; install it in the QGIS Python environment if needed. If `import scipy` fails with **NumPy ufunc** errors, upgrade both: `python -m pip install -U numpy scipy` using QGIS’s Python.

**Scripts folder:** Only install Processing tool scripts there — extra `.py` files in that folder are loaded as tools and can break the toolbox on import.

### Algorithm id
`tree_peaks_tiede_style`

### Workflow (summary)
1. **Preliminary mask:** discard cells below **Minimum canopy height** (and invalid/nodata).
2. **Per cell:** **CD = min(a + b·TH, max crown diameter cap)**, **r_search_m = CD/2** (m).
3. **Coarse candidates:** scipy `maximum_filter` on a metric disk of radius **max CD/2** (pixel footprint).
4. **Peaks:** among coarse candidates, cells that are local maxima in a disk of radius **r_search_m** at that cell; equal height neighbors disqualify unless this cell is **first in row-major order** `(row, col)` among ties in that disk.

### Inputs / outputs
- **Input raster (DSM or canopy height CHM/nDOM)** (RasterLayer)
- **DTM (optional)** — if set, **nDOM = DSM − DTM** (first input treated as DSM)
- **Extent (optional)** — limit processing to this rectangle (map/project CRS; must overlap the raster). Applied after nDOM when DTM is used.
- **Crown model: a** in CD = a + b·TH (meters) (Number, default 1.54)
- **Crown model: b** in CD = a + b·TH (Number, default 0.123)
- **Maximum crown diameter cap (meters)** (Number, default 10.0)
- **Minimum canopy height — preliminary mask (meters)** (Number, default 5.0)
- **Edge padding (meters)** (Number, default 10.0) — pad grid for full metric-disk sampling at borders
- **Output tree peaks** (VectorDestination, required) — attributes: `peak_id`, `tree_height`, `cd_eff_m`, `r_search_m`

### Parameters (defaults)
| Id | Type | Default |
|----|------|--------|
| `INPUT_RASTER` | RasterLayer | — |
| `INPUT_DTM` | RasterLayer | optional |
| `EXTENT` | Extent | optional |
| `COEFF_A` | Double | 1.54 |
| `COEFF_B` | Double | 0.123 |
| `MAX_CROWN_DIAMETER_M` | Double | 10.0 |
| `MIN_TREE_HEIGHT_M` | Double | 5.0 |
| `EDGE_PAD_M` | Double | 10.0 |
| `OUTPUT_PEAKS` | VectorDestination | — |

**Edge padding:** `EDGE_PAD_M` expands the grid on all sides by this distance (converted to pixels per axis). Padded cells use a height below the minimum canopy threshold so Tiede disks are fully sampled at the raster border; only the original extent is written to outputs. Set `0` to disable padding (legacy behaviour).

### Usage
- Python:
```python
import processing
processing.run('tree_peaks_tiede_style', {
    'INPUT_RASTER': '/path/to/dsm_or_chm.tif',
    'INPUT_DTM': None,  # or path to DTM for nDOM
    'EXTENT': None,  # or QgsRectangle / extent string from context
    'COEFF_A': 1.54,
    'COEFF_B': 0.123,
    'MAX_CROWN_DIAMETER_M': 10.0,
    'MIN_TREE_HEIGHT_M': 5.0,
    'EDGE_PAD_M': 10.0,
    'OUTPUT_PEAKS': '/path/to/peaks.gpkg',
})
```

If `processing.run` cannot find the id, use the id shown in Processing ▶ tool help (script tools may appear as `script:…`).

### Breaking change
Earlier versions exposed **`tree_crowns_tiede_style`** with crown raster and `OUTPUT_SEEDS` / `OUTPUT_TIPS`. Saved models and scripts must be updated to **`tree_peaks_tiede_style`** and **`OUTPUT_PEAKS`** (optional coarse-candidate layer was removed).
