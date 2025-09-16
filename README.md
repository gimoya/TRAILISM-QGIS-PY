# TRAILISM QGIS Processing Scripts

Lightweight processing algorithms for QGIS focused on trail design and geometry work. This repo currently contains three tools:

- TRAILISM: Find Curves — Detects curve apex points along linework
- TRAILISM: Round Corner by Min/Max Radius — Inserts circular fillets at polyline corners
- TRAILISM: Cross-Section Profiles — Samples DEM cross-sections, finds slope intersections, computes per-station cut/fill, and (optionally) exports 1:1 SVGs


## Requirements

- QGIS 3.x (tested with recent 3.2x–3.3x builds)
- Python bundled with QGIS
- Optional (for SVG export in Cross-Section Profiles):
  - matplotlib ≥ 3.5
  - numpy ≥ 1.22

Install optional packages into the QGIS Python environment if you want SVGs:

```bash
python -m pip install matplotlib numpy
```


## Installation

You already have these linked into your Processing toolbox — nothing to do. For a fresh setup:

Option A — User scripts folder (simple):
1) In QGIS: Settings → Options → Processing → “Scripts” → note the folder path
2) Copy the three .py files from this repo into that folder
3) Restart QGIS. The tools appear in the Processing Toolbox as “TRAILISM: …”.

Option B — As part of a plugin/provider (advanced):
- Place these algorithms inside a custom processing provider and load it as a QGIS plugin. This gives you versioning and distribution controls.


## Tools overview

1) TRAILISM: Find Curves
- Input: line layer
- Output: point layer marking curve apexes (attributes copied from input)
- Notes: works on vertex geometry; requires ≥ 3 points per feature

2) TRAILISM: Round Corner by Min/Max Radius
- Input: line layer
- Key parameters: min/max radius, max chord length, optional berm offset and threshold
- Output: new line layer with fillet arcs; optional berms line layer

3) TRAILISM: Cross-Section Profiles
- Inputs: trail centerline (line) and DEM raster
- Key parameters: station spacing, width (field or default), profile half-range R, side slope angle β1, DEM offset, axis dz or full cut/fill modes, SVG options
- Outputs:
  - Points: axis center, edges, and C points with per-station cut/fill areas
  - Lines: BC segments, chains for trail edges and C points, optional transect lines for plotted stations
  - Optional: 1:1 SVG cross-sections saved next to the QGIS project


## Using in QGIS (GUI)

- Open Processing Toolbox and search for “TRAILISM”.
- Each tool has inline help; parameters mirror the per-script docs in this repo.
- Use the Model Designer for batch runs (e.g., multiple layers or parameter sweeps).


## Using from the Python console

Algorithm ids depend on how you installed them. Use the toolbox tooltip or list algorithms to get the exact id, then run:

```python
import processing
processing.algorithmHelp('<algorithm_id>')
# Example (id will differ in your setup):
# processing.run('roundcornermaxradius', {'INPUT': layer, 'MIN_RADIUS': 3.0, 'MAX_RADIUS': 6.0, 'MAX_SEG_LEN': 0.2, 'DRAW_BERMS': False, 'OUTPUT': 'memory:'})
```


## Sidecar docs

See one-pagers next to each script for precise parameters, outputs, and tips:

- TRAILISM__Curve_Detection.md
- TRAILISM__Round_Corner_Linestrings_Radius.md
- TRAILISM__Transect_Profiles.md


## License

The Python sources contain GPL-2.0-or-later headers. This repository is distributed under the same terms.


## Notes and caveats

- Geometry is processed in the layer’s CRS. For DEM sampling, the tool transforms to the DEM’s CRS internally.
- SVG export is capped to 100 files per run to avoid accidental floods.
- “Find Curves” is experimental and depends on vertex density/quality; consider densifying/smoothing upstream for best results.

