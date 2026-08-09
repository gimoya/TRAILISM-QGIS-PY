## TRAILISM: Corner Fillets (min/max radius)

**Bisector fillets** (classical CAD-style fillets): each corner becomes a circular arc **inside** the turn. The arc center sits on the **interior angle bisector**; the arc is **tangent to both incoming and outgoing segments** (tangent points lie on those legs—not circles drawn on the vertex itself).

- **Cleanup:** You can enforce a minimum distance between vertices along the polyline so short noise segments are removed before filleting.
- **Radius:** The geometric limit per corner matches **TRAILISM: Round Corner by Min/Max Radius**: turning angle and leg lengths set how large a fillet can be. User Min/Max radius clamp that value. If the requested **minimum** radius is larger than geometry allows, the fillet uses **what the corner can support** (no overshoot past the segment).
- **Discretization:** Arcs are broken into straight segments with chord length ≤ **Max chord length**.
- **Berms:** Optional extra lines: offset arcs **outside** the fillet when radius is within your threshold (same pattern as the older round-corner berm output).

For **vertex-centered** circles and **outer common tangents**, use **TRAILISM: Corner Arcs (vertex circles, outside)**.

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
