## TRAILISM: Corner Fillets (min/max radius)

TRAILISM ALGORITHM (bisector fillets):\n\n- Optional: cleans vertices that are closer than a minimum spacing along the polyline.\n- At each corner, builds a **fillet**: the arc lies inside the turn; its center sits on the **interior angle bisector** (incenter-based), and the arc is **tangent to both legs** (tangent points on the adjacent segments).\n- Fillet radius uses the same geometric allowance as “Round corner by min/max radius”, then clamps to your Min/Max radius.\n- Arcs are approximated with straight segments so chord length ≤ Max chord length.\n- Optional berm output: offset arc outside the fillet (same as the original round-corner tool).\n\nThis is **not** vertex-centered rounding. For circles centered on each vertex with outer tangents, use “TRAILISM: Corner Arcs (vertex circles, outside)”.

### Algorithm id
`corner_fillets_min_max_radius`

### Inputs
- **Input layer** (FeatureSource)
- **Min vertex spacing before filleting (map units)** (Number)
- **Min radius (0 = no minimum)** (Number)
- **Max radius (0 = no maximum)** (Number)
- **Max chord length for arcs (0 = default tolerance)** (Number)
- **Draw berms as offset arc (for curves with radius ≤ threshold)** (Boolean)
- **Berm offset distance (map units)** (Number)
- **Berm radius threshold (draw berm when fillet radius ≤ value)** (Number)
- **Min berm length (map units); only output berms ≥ this length** (Number)

### Parameters
- **Input layer** (FeatureSource)
- **Min vertex spacing before filleting (map units)** (Number, default 1.0)
- **Min radius (0 = no minimum)** (Number, default 3.0)
- **Max radius (0 = no maximum)** (Number, default 6.0)
- **Max chord length for arcs (0 = default tolerance)** (Number, default 0.5)
- **Draw berms as offset arc (for curves with radius ≤ threshold)** (Boolean, default False)
- **Berm offset distance (map units)** (Number, default 1.0)
- **Berm radius threshold (draw berm when fillet radius ≤ value)** (Number, default 6.0)
- **Min berm length (map units); only output berms ≥ this length** (Number, default 5.0)
- **Berms layer (optional)** (FeatureSink)
- **Curved lines (filleted)** (FeatureSink)

### Outputs
- **Berms layer (optional)**
- **Curved lines (filleted)**

### Usage
- Python:
```python
import processing
processing.run('corner_fillets_min_max_radius', {
    'INPUT': ...,
    'MIN_VERTEX_SPACING': ...,
    'MIN_RADIUS': ...,
    'MAX_RADIUS': ...,
    'MAX_SEG_LEN': ...,
    'DRAW_BERMS': ...,
    'BERM_OFFSET': ...,
    'BERM_RADIUS_THRESHOLD': ...,
    'MIN_BERM_LENGTH': ...,
    'OUTPUT_BERMS': ...,
    'OUTPUT': ...,
    'OUTPUT': 'memory:'
})
```
