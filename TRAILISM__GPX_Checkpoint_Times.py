"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

import math
import os
import glob
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFile,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingException,
    QgsProcessingContext,
    QgsProcessingLayerPostProcessorInterface,
    QgsFeature,
    QgsFields,
    QgsField,
    QgsGeometry,
    QgsPoint,
    QgsLineString,
    QgsWkbTypes,
    QgsSingleSymbolRenderer,
    QgsLineSymbol,
    QgsProject,
)

GPX_NAMESPACES = (
    'http://www.topografix.com/GPX/1/1',
    'http://www.topografix.com/GPX/1/0',
)


def _haversine_m(lat1, lon1, lat2, lon2):
    """Distance in meters between two WGS84 points (lat/lon in degrees)."""
    R = 6371000
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


# Common GPX time formats (ISO 8601, RFC 3339, and variants from Garmin, Strava, etc.)
_TIME_FORMATS = [
    ('%Y-%m-%dT%H:%M:%S.%fZ', True),
    ('%Y-%m-%dT%H:%M:%S.%f+00:00', True),
    ('%Y-%m-%dT%H:%M:%SZ', True),
    ('%Y-%m-%dT%H:%M:%S+00:00', True),
    ('%Y-%m-%dT%H:%M:%S.%f', False),
    ('%Y-%m-%dT%H:%M:%S', False),
    ('%Y-%m-%d %H:%M:%S.%f', False),
    ('%Y-%m-%d %H:%M:%S', False),
    ('%Y%m%dT%H%M%SZ', True),
    ('%Y%m%dT%H%M%S', False),
    ('%Y%m%dT%H%M%S.%fZ', True),
    ('%d/%m/%Y %H:%M:%S', False),
    ('%d/%m/%Y %H:%M:%S.%f', False),
    ('%d.%m.%Y %H:%M:%S', False),
    ('%d-%m-%Y %H:%M:%S', False),
    ('%d-%m-%Y %H:%M:%S.%f', False),
    ('%m/%d/%Y %H:%M:%S', False),
    ('%m-%d-%Y %H:%M:%S', False),
]


def _parse_time(s):
    """Parse datetime string; return datetime or None. Tries common GPX formats."""
    s = (str(s) or '').strip()
    if not s:
        return None
    # Normalize Z and common offsets for fromisoformat
    s_norm = s.replace('Z', '+00:00').replace('+0000', '+00:00').replace('-0000', '-00:00')
    try:
        dt = datetime.fromisoformat(s_norm)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt
    except (ValueError, TypeError):
        pass
    for fmt, has_tz in _TIME_FORMATS:
        try:
            dt = datetime.strptime(s, fmt)
            if not has_tz and dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            return dt
        except (ValueError, TypeError):
            continue
    return None


def _duration_hhmmss(start_dt, end_dt):
    """Return duration as hh:mm:ss string."""
    if start_dt is None or end_dt is None:
        return None
    sec = int((end_dt - start_dt).total_seconds())
    if sec < 0:
        return None
    h, r = divmod(sec, 3600)
    m, s = divmod(r, 60)
    return f'{h:02d}:{m:02d}:{s:02d}'


def _point_wkt(lon, lat):
    """WKT for a point (lon, lat in degrees)."""
    return f'POINT ({lon} {lat})' if lon is not None and lat is not None else None


def _gpx_ns_from_root(root):
    """Return namespace map from root tag. Root must have namespace (caller validates)."""
    ns_uri = root.tag[1:].split('}', 1)[0]
    return {'gpx': ns_uri}


def _local_name(tag):
    """Return local name (without namespace) from Element tag."""
    if tag.startswith('{'):
        return tag.split('}', 1)[1]
    return tag


def _get_time_text(trkpt, ns):
    """Get time string from trkpt. Tries namespaced find, then local name match."""
    time_el = trkpt.find('gpx:time', ns)
    if time_el is not None and time_el.text:
        return (time_el.text or '').strip()
    for child in trkpt:
        if _local_name(child.tag) == 'time' and child.text:
            return child.text.strip()
    return ''


class _InvalidGPXError(Exception):
    """Raised when a GPX file does not meet requirements (e.g. missing namespace)."""
    pass


def _parse_gpx_tracks(path):
    """
    Return list of tracks; each track is list of (lat, lon, time_dt).
    Requires proper GPX namespace (GPX 1.0 or 1.1). Raises _InvalidGPXError if
    the root element has no xmlns (root tag does not start with '{').
    """
    path_str = os.path.abspath(os.path.normpath(str(path).strip()))
    if not path_str or not os.path.isfile(path_str):
        return []

    tree = ET.parse(path_str)
    root = tree.getroot()
    if root is None:
        return []

    # Proper GPX must declare a namespace; root tag should be {uri}gpx
    if not root.tag.startswith('{'):
        raise _InvalidGPXError(
            os.path.basename(path_str),
            "GPX root element has no XML namespace (xmlns). "
            "GPX 1.0/1.1 must declare e.g. xmlns=\"http://www.topografix.com/GPX/1/1\" on the root <gpx> element."
        )

    def collect_pt(el, namespace):
        lat = el.get('lat')
        lon = el.get('lon')
        if lat is None or lon is None:
            return None
        try:
            lat_f = float(lat)
            lon_f = float(lon)
        except ValueError:
            return None
        time_str = _get_time_text(el, namespace)
        time_dt = _parse_time(time_str)
        return (lat_f, lon_f, time_dt)

    tracks = []
    root_ns = _gpx_ns_from_root(root)
    ns_candidates = [root_ns] + [{'gpx': uri} for uri in GPX_NAMESPACES if uri != root_ns.get('gpx')]
    for ns in ns_candidates:
        for trk in root.findall('.//gpx:trk', ns):
            pts = []
            for trkseg in trk.findall('gpx:trkseg', ns):
                for trkpt in trkseg.findall('gpx:trkpt', ns):
                    pt = collect_pt(trkpt, ns)
                    if pt is not None:
                        pts.append(pt)
            if pts:
                tracks.append(pts)
        for rte in root.findall('.//gpx:rte', ns):
            pts = []
            for rtept in rte.findall('gpx:rtept', ns):
                pt = collect_pt(rtept, ns)
                if pt is not None:
                    pts.append(pt)
            if pts:
                tracks.append(pts)
        if tracks:
            break
    return tracks


def _build_edges(track_pts):
    """
    Build edge list from track points. Each edge: (start, end, dist_m, time_s).
    start/end = (lat, lon); dist_m = haversine; time_s = duration to next point.
    """
    n = len(track_pts)
    edges = []
    for i in range(n - 1):
        lat0, lon0, t0 = track_pts[i]
        lat1, lon1, t1 = track_pts[i + 1]
        dist = _haversine_m(lat0, lon0, lat1, lon1)
        time_s = (t1 - t0).total_seconds() if t0 and t1 else 0.0
        edges.append(((lat0, lon0), (lat1, lon1), dist, time_s))
    return edges


def _find_all_snaps(track_pts, ref_pts, max_dist_m):
    """
    For each ref point, find the track point closest to it within max_dist_m.
    Returns [(idx, lat, lon, time_dt), ...] sorted by idx, one per ref that was matched.
    """
    if not track_pts or not ref_pts:
        return []
    result = []
    used_track = set()
    for ref_lat, ref_lon in ref_pts:
        best_idx = None
        best_dist = float('inf')
        for i, (lat, lon, t) in enumerate(track_pts):
            if i in used_track:
                continue
            d = _haversine_m(lat, lon, ref_lat, ref_lon)
            if d <= max_dist_m and d < best_dist:
                best_dist = d
                best_idx = i
        if best_idx is not None:
            used_track.add(best_idx)
            lat, lon, t = track_pts[best_idx]
            result.append((best_idx, lat, lon, t))
    return sorted(result, key=lambda x: x[0])


def _pair_segments(start_snaps, end_snaps):
    """
    Pair by position: start1-end1, start2-end2, ... Only valid if end_i > start_i.
    Returns [(start_idx, end_idx), ...] and counts (n_start_orphan, n_end_orphan, n_invalid).
    """
    n = min(len(start_snaps), len(end_snaps))
    pairs = []
    n_invalid = 0
    for i in range(n):
        s_idx = start_snaps[i][0]
        e_idx = end_snaps[i][0]
        if e_idx > s_idx:
            pairs.append((s_idx, e_idx))
        else:
            n_invalid += 1
    n_start_orphan = len(start_snaps) - n
    n_end_orphan = len(end_snaps) - n
    return pairs, n_start_orphan, n_end_orphan, n_invalid


def _segment_from_edges(edges, s_idx, e_idx):
    """
    From edges, compute segment [s_idx, e_idx]: dist, time, and pts [(lat,lon),...].
    Uses edges s_idx to e_idx-1.
    """
    dist = sum(edges[i][2] for i in range(s_idx, e_idx))
    time_s = sum(edges[i][3] for i in range(s_idx, e_idx))
    pts = [edges[s_idx][0]] + [edges[i][1] for i in range(s_idx, e_idx)]
    return dist, time_s, pts


def _get_output_layer(context, dest_id):
    """Resolve output layer from context (temp store or project)."""
    return context.temporaryLayerStore().mapLayer(dest_id) or (
        QgsProject.instance().mapLayer(dest_id) if context.project() else None
    )


def _style_line_layer(layer, color_hex, width):
    """Apply simple line symbol and trigger repaint."""
    if layer is None:
        return
    sym = QgsLineSymbol.createSimple({"color": color_hex, "width": str(width)})
    layer.setRenderer(QgsSingleSymbolRenderer(sym))
    layer.triggerRepaint()


class _RedStylePostProcessor(QgsProcessingLayerPostProcessorInterface):
    def postProcessLayer(self, layer, context, feedback):
        _style_line_layer(layer, "#dc3232", 1.8)


class _GreenStylePostProcessor(QgsProcessingLayerPostProcessorInterface):
    def postProcessLayer(self, layer, context, feedback):
        _style_line_layer(layer, "#00c800", 1.2)


class GpxCheckpointTimesAlgorithm(QgsProcessingAlgorithm):
    """
    Start/end point layers (WGS84), folder of GPX files. Pairs each start with next end along track.
    Output: one feature per matched segment + original tracks. Edge-based: dist/time from summed edges.
    """

    INPUT_START = 'INPUT_START'
    INPUT_END = 'INPUT_END'
    INPUT_GPX_FILES = 'INPUT_GPX_FILES'
    MAX_DISTANCE_M = 'MAX_DISTANCE_M'
    OUTPUT_LINES = 'OUTPUT_LINES'
    OUTPUT_ORIGINAL_TRACKS = 'OUTPUT_ORIGINAL_TRACKS'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return GpxCheckpointTimesAlgorithm()

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def name(self):
        return 'gpx_checkpoint_times'

    def displayName(self):
        return self.tr('TRAILISM: GPX Checkpoint Times')

    def shortHelpString(self):
        return self.tr(
            'Start/end point layers (WGS84), folder of GPX files. Pairs each start with next end along track.\n'
            'Output: matched segments (one per start-end pair) + original tracks for comparison. '
            'Attributes: match_info, start_wkt, end_wkt, segment_dist, start/end times, duration.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_START,
                self.tr('Start point layer (x/y in WGS84)'),
                [QgsProcessing.TypeVectorPoint]
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_END,
                self.tr('End point layer (x/y in WGS84)'),
                [QgsProcessing.TypeVectorPoint]
            )
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_GPX_FILES,
                self.tr('Folder containing GPX files'),
                behavior=QgsProcessingParameterFile.Folder,
                fileFilter='GPX files (*.gpx);;All files (*.*)'
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_DISTANCE_M,
                self.tr('Max distance (m) for track points to start/end'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=15.0,
                minValue=0.0,
                optional=True
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_ORIGINAL_TRACKS,
                self.tr('Output original GPX tracks (red)')
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_LINES,
                self.tr('Output matched segments (green)')
            )
        )

    def _crs_is_wgs84(self, crs):
        return crs is not None and crs.authid() == 'EPSG:4326'

    def processAlgorithm(self, parameters, context, feedback):
        start_source = self.parameterAsSource(parameters, self.INPUT_START, context)
        end_source = self.parameterAsSource(parameters, self.INPUT_END, context)
        folder = self.parameterAsString(parameters, self.INPUT_GPX_FILES, context).strip()
        if not folder or not os.path.isdir(folder):
            raise QgsProcessingException(self.tr('Invalid or empty folder.'))
        gpx_paths = [p for p in glob.glob(os.path.join(folder, '*.gpx')) if os.path.isfile(p)]
        if start_source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_START))
        if end_source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_END))
        if not gpx_paths:
            raise QgsProcessingException(self.tr('No GPX files found in folder.'))

        max_dist_m = max(0.0, self.parameterAsDouble(parameters, self.MAX_DISTANCE_M, context))

        for name, src in [('Start', start_source), ('End', end_source)]:
            if not self._crs_is_wgs84(src.sourceCrs()):
                raise QgsProcessingException(
                    self.tr('{} layer must be WGS 84 (EPSG:4326). Current: {}').format(
                        name, src.sourceCrs().authid() or 'unknown'
                    )
                )

        start_pts = [(f.geometry().asPoint().y(), f.geometry().asPoint().x())
                    for f in start_source.getFeatures()
                    if f.geometry() and not f.geometry().isEmpty()]
        end_pts = [(f.geometry().asPoint().y(), f.geometry().asPoint().x())
                   for f in end_source.getFeatures()
                   if f.geometry() and not f.geometry().isEmpty()]

        if not start_pts:
            raise QgsProcessingException(self.tr('Start point layer has no valid points.'))
        if not end_pts:
            raise QgsProcessingException(self.tr('End point layer has no valid points.'))

        line_fields = QgsFields()
        line_fields.append(QgsField('source_file', QVariant.String))
        line_fields.append(QgsField('match_info', QVariant.String))
        line_fields.append(QgsField('start_wkt', QVariant.String))
        line_fields.append(QgsField('end_wkt', QVariant.String))
        line_fields.append(QgsField('segment_dist', QVariant.Double))
        line_fields.append(QgsField('start_time', QVariant.String))
        line_fields.append(QgsField('end_time', QVariant.String))
        line_fields.append(QgsField('duration', QVariant.String))

        (line_sink, line_dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT_LINES, context,
            line_fields, QgsWkbTypes.LineString, start_source.sourceCrs()
        )
        if line_sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT_LINES))

        orig_fields = QgsFields()
        orig_fields.append(QgsField('source_file', QVariant.String))
        orig_fields.append(QgsField('match_info', QVariant.String))
        (orig_sink, orig_dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT_ORIGINAL_TRACKS, context,
            orig_fields, QgsWkbTypes.LineString, start_source.sourceCrs()
        )
        if orig_sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT_ORIGINAL_TRACKS))

        gpx_tracks = []
        for p in gpx_paths:
            try:
                tracks = _parse_gpx_tracks(p)
                gpx_tracks.append((p, tracks))
            except _InvalidGPXError as e:
                raise QgsProcessingException(
                    self.tr('Invalid GPX "{}": {}').format(e.args[0], e.args[1])
                )
            except Exception as e:
                raise QgsProcessingException(
                    self.tr('Failed to parse GPX "{}": {}').format(os.path.basename(p), str(e))
                )
        total_tracks = sum(len(t) for _, t in gpx_tracks)
        done = 0

        for gpx_path, tracks in gpx_tracks:
            if feedback.isCanceled():
                break
            basename = os.path.basename(gpx_path)
            for track_pts in tracks:
                if len(track_pts) < 2:
                    continue

                edges = _build_edges(track_pts)
                start_snaps = _find_all_snaps(track_pts, start_pts, max_dist_m)
                end_snaps = _find_all_snaps(track_pts, end_pts, max_dist_m)
                pairs, n_start_orphan, n_end_orphan, n_invalid = _pair_segments(start_snaps, end_snaps)

                n_ok, n_start, n_end = len(pairs), len(start_snaps), len(end_snaps)
                if n_ok == 0 and n_start == 0 and n_end == 0:
                    track_match = "FAIL no_snaps"
                elif n_ok == 0 and n_start == 0:
                    track_match = "FAIL no_start"
                elif n_ok == 0 and n_end == 0:
                    track_match = "FAIL no_end"
                else:
                    parts = [f"{n_ok} OK" if n_ok else "0"]
                    if n_start_orphan:
                        parts.append(f"{n_start_orphan} no_end")
                    if n_end_orphan:
                        parts.append(f"{n_end_orphan} no_start")
                    if n_invalid:
                        parts.append(f"{n_invalid} invalid")
                    track_match = " ".join(parts)

                for s_idx, e_idx in pairs:
                    dist, time_s, pts = _segment_from_edges(edges, s_idx, e_idx)
                    start_lat, start_lon = pts[0]
                    end_lat, end_lon = pts[-1]
                    start_time = track_pts[s_idx][2]
                    end_time = track_pts[e_idx][2]
                    duration_str = _duration_hhmmss(start_time, end_time)
                    match_info = "OK {}→{} ({:.0f}m {:.0f}s)".format(s_idx, e_idx, dist, time_s)

                    qpts = [QgsPoint(lon, lat) for lat, lon in pts]
                    feat = QgsFeature(line_fields)
                    feat.setGeometry(QgsGeometry(QgsLineString(qpts)))
                    feat.setAttributes([
                        basename, match_info,
                        _point_wkt(start_lon, start_lat), _point_wkt(end_lon, end_lat),
                        dist,
                        start_time.isoformat() if start_time else None,
                        end_time.isoformat() if end_time else None,
                        duration_str
                    ])
                    line_sink.addFeature(feat)

                # Always add orig track to output, regardless of match status (even when 0 pairs)
                orig_pts = [QgsPoint(pt[1], pt[0]) for pt in track_pts]
                orig_feat = QgsFeature(orig_fields)
                orig_feat.setGeometry(QgsGeometry(QgsLineString(orig_pts)))
                orig_feat.setAttributes([basename, track_match])
                orig_sink.addFeature(orig_feat)

                done += 1
                if total_tracks > 0:
                    feedback.setProgress(int(100 * done / total_tracks))

        # Pre-style immediately if layer in context; post-processor handles late load (e.g. file output)
        _style_line_layer(_get_output_layer(context, orig_dest_id), "#dc3232", 1.8)
        _style_line_layer(_get_output_layer(context, line_dest_id), "#00c800", 1.2)
        if context.project():
            for dest_id, pp, name, output_key, sort_key in [
                (orig_dest_id, _RedStylePostProcessor(), self.tr("Original GPX tracks"), self.OUTPUT_ORIGINAL_TRACKS, 0),
                (line_dest_id, _GreenStylePostProcessor(), self.tr("Matched segments"), self.OUTPUT_LINES, 1),
            ]:
                details = QgsProcessingContext.LayerDetails(name, context.project(), output_key)
                details.layerSortKey = sort_key
                details.setPostProcessor(pp)
                context.addLayerToLoadOnCompletion(dest_id, details)

        return {self.OUTPUT_ORIGINAL_TRACKS: orig_dest_id, self.OUTPUT_LINES: line_dest_id}
