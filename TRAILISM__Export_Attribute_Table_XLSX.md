## TRAILISM: Export Attribute Table to XLSX

TRAILISM ALGORITHM:\n\nExports the attribute table of a vector layer to Excel (.xlsx).\nOnly data columns — no geometry, no WKT, no layer metadata.\n\nMinimal input:\n- Layer (defaults to the active layer when available)\n- Output folder: last-used folder first, else layer source folder\n\nFile name is always the layer name + .xlsx.

### Algorithm id
`export_attribute_table_xlsx`

### Inputs
- **Input layer** (VectorLayer)
- **Output folder** (File)

### Parameters
- **Input layer** (VectorLayer)
- **Output folder** (File, default default_dir)

### Usage
- Python:
```python
import processing
processing.run('export_attribute_table_xlsx', {
    'INPUT': ...,
    'OUTPUT_DIR': ...,
    'OUTPUT': 'memory:'
})
```
