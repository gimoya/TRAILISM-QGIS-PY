## TRAILISM: Round Corner with Fixed Radius

Rounds corners with a fixed radius applied to all vertices. The radius is applied to all corners, trimmed to fit within adjacent segments.

### Algorithm id
`roundcornerfixedradius` (check your Processing id via tooltip)

### Inputs
- **INPUT**: Vector line layer

### Parameters
- **RADIUS** (Double, default 5.0): Fixed radius applied to all corners (map units). If the radius is too large for a corner, it will be trimmed to fit within the adjacent segments.
- **MAX_SEG_LEN** (Double, default 0.2): Maximum chord length for arc linearization (map units). 0 uses a deviation tolerance of 0.1.

### Outputs
- **OUTPUT**: Curved linestrings with fixed radius fillets applied.

### Behavior
- Fixed radius applied to every corner vertex.
- Per-vertex radius cap: r_allowed = 0.5 · min(Lprev, Lnext) · tan(θ/2). Final radius is min(fixed_radius, r_allowed).
- Arc center picked on the angle bisector. The shorter sweep matching φ = π − θ is used.
- Linearization uses `MAX_SEG_LEN` or deviation 0.1 map units when 0.
- Output layer is styled white (0.3 mm) by default.

### Usage
- GUI: Processing Toolbox → "TRAILISM: Round Corner with Fixed Radius".
- Python:
```python
import processing
processing.run('roundcornerfixedradius', {
    'INPUT': input_layer,
    'RADIUS': 5.0,
    'MAX_SEG_LEN': 0.2,
    'OUTPUT': 'memory:'
})
```

### Tips
- Use a smaller fixed radius for tighter, more consistent curves.
- Increase MAX_SEG_LEN for smoother arcs with fewer segments.
- The algorithm automatically trims radius when geometry constraints require it.

