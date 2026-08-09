## TRAILISM: GPX Checkpoint Times

Start/end point layers (WGS84), folder of GPX files. Pairs each start with next end along track.\nOutput: matched segments (one per start-end pair) + original tracks for comparison. Attributes: match_info, start_wkt, end_wkt, segment_dist, start/end times, duration.

### Algorithm id
`gpx_checkpoint_times`

### Inputs
- **Start point layer (x/y in WGS84)** (FeatureSource)
- **End point layer (x/y in WGS84)** (FeatureSource)
- **Folder containing GPX files** (File)
- **Max distance (m) for track points to start/end** (Number)

### Parameters
- **Start point layer (x/y in WGS84)** (FeatureSource)
- **End point layer (x/y in WGS84)** (FeatureSource)
- **Folder containing GPX files** (File)
- **Max distance (m) for track points to start/end** (Number, default 15.0)
- **Output original GPX tracks (red)** (FeatureSink)
- **Output matched segments (green)** (FeatureSink)

### Outputs
- **Output original GPX tracks (red)**
- **Output matched segments (green)**

### Usage
- Python:
```python
import processing
processing.run('gpx_checkpoint_times', {
    'INPUT_START': ...,
    'INPUT_END': ...,
    'INPUT_GPX_FILES': ...,
    'MAX_DISTANCE_M': ...,
    'OUTPUT_ORIGINAL_TRACKS': ...,
    'OUTPUT_LINES': ...,
    'OUTPUT': 'memory:'
})
```
