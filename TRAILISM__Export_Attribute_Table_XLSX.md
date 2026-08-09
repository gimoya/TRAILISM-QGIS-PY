## TRAILISM: Export Attribute Table to XLSX

One-shot export of a vector layer attribute table to Excel (`.xlsx`). Data columns only — no geometry, no WKT, no metadata. File name = layer name.

### Algorithm id
`export_attribute_table_xlsx`

### Inputs
- **Input layer** (VectorLayer)
- **Output folder** (Folder, optional; default = last-used, else layer source dir)

### Outputs
- `{layer_name}.xlsx` in the chosen folder

### Usage
Processing Toolbox → TRAILISM → Export Attribute Table to XLSX

```python
import processing
processing.run('script:export_attribute_table_xlsx', {
    'INPUT': 'my_layer_id',
    'OUTPUT_DIR': r'C:\exports',
})
```
