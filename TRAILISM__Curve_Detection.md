## TRAILISM: Find Curves

### What it does
Detects curve sections in a line layer and marks likely apex points based on arc-length to chord-length patterns between turn points.

### Algorithm id
`find curves` (exact id depends on how you registered the provider)

### Inputs
- **INPUT**: Vector line layer

### Outputs
- **OUTPUT**: Point features at detected curve apex locations. Attributes are copied from the input feature.

### Parameters and behavior
- Works on polyline vertices. Needs ≥ 3 vertices to analyze a section.
- Internally identifies turn points via direction change, then analyzes sub-sections using arc/chord ratio trends to find apex.

### Usage
- GUI: Processing Toolbox → “TRAILISM: Find Curves”.
- Python console (id may differ in your setup):
```python
import processing
processing.run('find curves', {
    'INPUT': input_layer,
    'OUTPUT': 'memory:'
})
```

### Tips
- Densify input lines if vertices are sparse; apex detection benefits from consistent vertex spacing.
- Remove noisy spikes before running to avoid false positives.

### Notes
- Output type should be points. If your environment errors about geometry type, ensure the tool is registered so the sink is a point layer.


