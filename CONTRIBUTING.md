# Contributing

## Adding a new QGIS Processing tool

1) Add your `.py` to the repo and register it in QGIS as usual.
2) Create docs with the scaffold (generates a sidecar and a docs page, and updates nav):
```bash
python scripts/scaffold_tool.py --file TRAILISM__Your_Tool.py --validate
```
3) Push to `main`. GitHub Actions builds and deploys the site.

Notes:
- The docs site includes sidecars directly. Keep the sidecar next to the script.
- If you rename a script, rename its sidecar and update the docs page slug (or re-run the scaffold with `--slug`).
- Local preview: `pip install -r requirements-docs.txt && mkdocs serve`.

## Style
- Prefer concise, task-focused docs: What it does, Inputs, Parameters, Outputs, Usage, Tips.
- Keep examples simple; prioritize readability over performance.
