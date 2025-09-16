## TRAILISM: Round Corner by Min/Max Radius

Smooths polyline corners by inserting a single circular arc per vertex inside the smaller angle. Optionally draws outward offset “berms” for small-radius fillets.

### Algorithm id
`roundcornermaxradius` (check your Processing id via tooltip)

### Inputs
- **INPUT**: Vector line layer

### Parameters
- **MIN_RADIUS** (Double, default 3.0): Minimum fillet radius. If local geometry cannot support it, tangent points are clamped to segment endpoints to enforce the minimum.
- **MAX_RADIUS** (Double, default 6.0): Maximum fillet radius. 0 disables maximum.
- **MAX_SEG_LEN** (Double, default 0.2): Maximum chord length for arc linearization (map units). 0 uses a deviation tolerance of 0.1.
- **DRAW_BERMS** (Boolean, default False): If true, emit offset arcs for fillets whose radius ≤ threshold.
- **BERM_OFFSET** (Double, default 1.0): Offset distance for berm arcs.
- **BERM_RADIUS_THRESHOLD** (Double, default 6.0): Draw berm only when fillet radius ≤ this value.

### Outputs
- **OUTPUT**: Curved linestrings with fillets applied.
- **OUTPUT_BERMS** (optional): Separate layer with berm arcs when enabled.

### Behavior
- Per-vertex radius cap: r_allowed = 0.5 · min(Lprev, Lnext) · tan(θ/2). Final radius is clamped by min/max parameters.
- Arc center picked on the angle bisector. The shorter sweep matching φ = π − θ is used.
- Linearization uses `MAX_SEG_LEN` or deviation 0.1 map units when 0.
- Output layers are styled white (0.3 mm) by default.

### Usage
- GUI: Processing Toolbox → “TRAILISM: Round Corner by Min/Max Radius”.
- Python:
```python
import processing
processing.run('roundcornermaxradius', {
    'INPUT': input_layer,
    'MIN_RADIUS': 3.0,
    'MAX_RADIUS': 6.0,
    'MAX_SEG_LEN': 0.2,
    'DRAW_BERMS': False,
    'BERM_OFFSET': 1.0,
    'BERM_RADIUS_THRESHOLD': 6.0,
    'OUTPUT': 'memory:',
    'OUTPUT_BERMS': 'memory:'
})
```

### Tips
- If you see “lasso-like” artifacts, reduce MAX_RADIUS or increase MIN_RADIUS to a feasible range.
- Densify upstream if long segments cause overly large `r_allowed` compared to your intent.


