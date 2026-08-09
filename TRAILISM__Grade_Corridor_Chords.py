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

from __future__ import annotations

import math

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsLineSymbol,
    QgsPointXY,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingLayerPostProcessorInterface,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterPoint,
    QgsProcessingParameterRasterLayer,
    QgsRuleBasedRenderer,
    QgsUnitTypes,
    QgsWkbTypes,
)


def _get_output_layer(context, dest_id):
    return context.temporaryLayerStore().mapLayer(dest_id) or (
        context.project().mapLayer(dest_id) if context.project() else None
    )


def _style_path_tries_layer(layer):
    """Prominent: end reached / longest; muted: other tries."""
    if layer is None:
        return
    root = QgsRuleBasedRenderer.Rule(None)

    sym_end = QgsLineSymbol.createSimple({
        'color': '232,93,4,255',
        'width': '1.6',
        'capstyle': 'round',
        'joinstyle': 'round',
    })
    rule_end = QgsRuleBasedRenderer.Rule(sym_end)
    rule_end.setFilterExpression('"reaches_end" = 1')
    rule_end.setLabel('End reached')
    root.appendChild(rule_end)

    sym_long = QgsLineSymbol.createSimple({
        'color': '244,140,6,230',
        'width': '1.2',
        'capstyle': 'round',
        'joinstyle': 'round',
    })
    rule_long = QgsRuleBasedRenderer.Rule(sym_long)
    rule_long.setFilterExpression(
        '"length_m" = maximum("length_m") AND "reaches_end" <> 1'
    )
    rule_long.setLabel('Longest try')
    root.appendChild(rule_long)

    sym_other = QgsLineSymbol.createSimple({
        'color': '120,120,120,110',
        'width': '0.35',
        'capstyle': 'round',
        'joinstyle': 'round',
    })
    rule_other = QgsRuleBasedRenderer.Rule(sym_other)
    rule_other.setIsElse(True)
    rule_other.setLabel('Other tries')
    root.appendChild(rule_other)

    layer.setRenderer(QgsRuleBasedRenderer(root))
    layer.triggerRepaint()


class _PathTriesStylePostProcessor(QgsProcessingLayerPostProcessorInterface):
    def postProcessLayer(self, layer, context, feedback):
        _style_path_tries_layer(layer)


class GradeCorridorChordsAlgorithm(QgsProcessingAlgorithm):
    """
    Circle steepest-down orientation → L/R semicircle grade pick → edge flip.
    On dead-end: drop tip, rewind to interior node, try opposite exit side once.
    """

    DEM = 'DEM'
    CORRIDOR = 'CORRIDOR'
    START = 'START'
    END = 'END'
    END_BUFFER = 'END_BUFFER'
    STEP_LENGTH = 'STEP_LENGTH'
    SIDE = 'SIDE'
    GRADE_PCT = 'GRADE_PCT'
    OUTPUT = 'OUTPUT'
    LINE_OUTPUT = 'LINE_OUTPUT'
    FALLLINE_OUTPUT = 'FALLLINE_OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return GradeCorridorChordsAlgorithm()

    def name(self):
        return 'grade_corridor_chords'

    def displayName(self):
        return self.tr('TRAILISM: Grade Corridor Chords')

    def group(self):
        return self.tr('TRAILISM')

    def groupId(self):
        return 'trailism'

    def shortHelpString(self):
        return self.tr(
            'TRAILISM ALGORITHM:\n\n'
            'Walk downhill inside a corridor:\n'
            '1. Full circle at chord length → steepest-down orientation '
            '(ties → angular mid).\n'
            '2. UI L/R semicircle (0…180°) → closest |grade − target|.\n'
            '3. Ray inside → next node; ray exits → edge intersection, flip side, continue.\n'
            '4. Two consecutive short chords / no candidate → drop tip, rewind, '
            'opposite exit side (once per node/side).\n'
            '5. End when within UI buffer of end point, or rewind exhausts back at start.\n\n'
            'Output: node points + all try paths as LineStrings + debug '
            'fall-line ray and 180° pick fan (best path).'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.DEM, self.tr('DEM raster')))
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.CORRIDOR,
                self.tr('Corridor polygon'),
                [QgsProcessing.TypeVectorPolygon],
            )
        )
        self.addParameter(
            QgsProcessingParameterPoint(self.START, self.tr('Start point (click on map)'))
        )
        self.addParameter(
            QgsProcessingParameterPoint(self.END, self.tr('End point (click on map)'))
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.END_BUFFER,
                self.tr('End reach buffer (m)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=5.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.STEP_LENGTH,
                self.tr('Chord length (m)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=30.0,
                minValue=1.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterEnum(
                self.SIDE,
                self.tr('Side of steepest-down orientation'),
                options=[self.tr('Left'), self.tr('Right')],
                defaultValue=0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.GRADE_PCT,
                self.tr('Desired grade (%)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=8.0,
                minValue=0.001,
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Nodes'),
                type=QgsProcessing.TypeVectorPoint,
                optional=True,
                createByDefault=False,
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.LINE_OUTPUT,
                self.tr('Path tries (all LineStrings)'),
                type=QgsProcessing.TypeVectorLine,
                optional=True,
                createByDefault=False,
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.FALLLINE_OUTPUT,
                self.tr('Debug: fall-line + 180° pick fan'),
                type=QgsProcessing.TypeVectorLine,
                optional=True,
                createByDefault=False,
            )
        )

    @staticmethod
    def _point_at(origin: QgsPointXY, bearing_deg: float, dist: float) -> QgsPointXY:
        rad = math.radians(bearing_deg)
        return QgsPointXY(
            origin.x() + math.sin(rad) * dist,
            origin.y() + math.cos(rad) * dist,
        )

    @staticmethod
    def _segment_inside(corridor: QgsGeometry, a: QgsPointXY, b: QgsPointXY, eps: float = 1e-3) -> bool:
        seg = QgsGeometry.fromPolylineXY([a, b])
        if seg.isEmpty():
            return False
        outside = seg.difference(corridor)
        if outside is None or outside.isEmpty():
            return True
        return outside.length() < eps

    @staticmethod
    def _angular_midpoint_deg(bearings: list[float]) -> float:
        if not bearings:
            raise ValueError('empty bearings')
        sx = sum(math.sin(math.radians(b)) for b in bearings)
        sy = sum(math.cos(math.radians(b)) for b in bearings)
        if abs(sx) < 1e-15 and abs(sy) < 1e-15:
            return bearings[0] % 360.0
        return (math.degrees(math.atan2(sx, sy)) + 360.0) % 360.0

    def processAlgorithm(self, parameters, context, feedback):
        dem = self.parameterAsRasterLayer(parameters, self.DEM, context)
        if dem is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.DEM))

        corridor_source = self.parameterAsSource(parameters, self.CORRIDOR, context)
        if corridor_source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.CORRIDOR))

        chord_length = self.parameterAsDouble(parameters, self.STEP_LENGTH, context)
        end_buffer = self.parameterAsDouble(parameters, self.END_BUFFER, context)
        side = self.parameterAsEnum(parameters, self.SIDE, context)
        grade_pct = self.parameterAsDouble(parameters, self.GRADE_PCT, context)
        g_target = grade_pct / 100.0
        side_idx = side
        side_label = 'L' if side_idx == 0 else 'R'

        corridor_crs = corridor_source.sourceCrs()
        if corridor_crs.isGeographic() or corridor_crs.mapUnits() != QgsUnitTypes.DistanceMeters:
            raise QgsProcessingException('Corridor CRS must be a projected CRS in metres.')

        dem_crs = dem.crs()
        if dem_crs.isGeographic() or dem_crs.mapUnits() != QgsUnitTypes.DistanceMeters:
            raise QgsProcessingException('DEM CRS must be a projected CRS in metres.')

        dem_provider = dem.dataProvider()
        if dem_provider is None:
            raise QgsProcessingException('DEM has no data provider.')

        transform_context = context.transformContext()
        xform_to_dem = None
        if corridor_crs != dem_crs:
            xform_to_dem = QgsCoordinateTransform(corridor_crs, dem_crs, transform_context)

        cell_dem = max(abs(dem.rasterUnitsPerPixelX()), abs(dem.rasterUnitsPerPixelY()))
        if cell_dem <= 0:
            raise QgsProcessingException('Invalid DEM cell size.')
        if xform_to_dem is None:
            cell = cell_dem
        else:
            dem_extent = dem.extent()
            p0 = QgsPointXY(dem_extent.xMinimum(), dem_extent.yMaximum())
            p1 = QgsPointXY(dem_extent.xMinimum() + cell_dem, dem_extent.yMaximum())
            xform_from_dem = QgsCoordinateTransform(dem_crs, corridor_crs, transform_context)
            cell = max(xform_from_dem.transform(p0).distance(xform_from_dem.transform(p1)), 1e-6)

        arc_step_deg = math.degrees(cell / chord_length)
        if arc_step_deg <= 0:
            raise QgsProcessingException('Invalid arc step from DEM resolution.')
        feedback.pushInfo(
            f'Circle/semicircle step {arc_step_deg:.4f}° '
            f'(cell {cell:.3f} m, chord {chord_length:.1f} m), side={side_label}'
        )

        corridor_parts = []
        for feat in corridor_source.getFeatures():
            geom = feat.geometry()
            if geom is None or geom.isEmpty():
                continue
            corridor_parts.append(QgsGeometry(geom))
        if not corridor_parts:
            raise QgsProcessingException('No valid corridor polygon features.')
        corridor = corridor_parts[0]
        for part in corridor_parts[1:]:
            corridor = corridor.combine(part)
        if corridor is None or corridor.isEmpty():
            raise QgsProcessingException('Corridor union is empty.')

        corridor_boundary = QgsGeometry(corridor.constGet().boundary())

        def fail(msg: str):
            feedback.reportError(msg)
            raise QgsProcessingException(msg)

        def point_in_corridor(xy: QgsPointXY) -> bool:
            geom = QgsGeometry.fromPointXY(xy)
            return bool(corridor.contains(geom) or corridor.intersects(geom))

        start_xy = self.parameterAsPoint(parameters, self.START, context, corridor_crs)
        end_xy = self.parameterAsPoint(parameters, self.END, context, corridor_crs)
        if not point_in_corridor(start_xy):
            fail(
                f'Start point is outside the corridor '
                f'(x={start_xy.x():.3f}, y={start_xy.y():.3f}).'
            )
        if not point_in_corridor(end_xy):
            fail(
                f'End point is outside the corridor '
                f'(x={end_xy.x():.3f}, y={end_xy.y():.3f}).'
            )
        if start_xy.distance(end_xy) < 1e-6:
            fail('Start and end points are coincident.')

        nodata = None
        if dem_provider.sourceHasNoDataValue(1):
            nodata = dem_provider.sourceNoDataValue(1)

        def sample_z(xy: QgsPointXY):
            dem_xy = xy if xform_to_dem is None else xform_to_dem.transform(xy)
            val, ok = dem_provider.sample(dem_xy, 1)
            bad = (
                (not ok)
                or val is None
                or (isinstance(val, float) and (math.isnan(val) or math.isinf(val)))
                or (nodata is not None and val == nodata)
            )
            if bad:
                return None
            return float(val)

        def edge_hit(origin: QgsPointXY, bearing_deg: float, d_cap: float):
            full = self._point_at(origin, bearing_deg, d_cap)
            if self._segment_inside(corridor, origin, full):
                return None
            inter = QgsGeometry.fromPolylineXY([origin, full]).intersection(corridor_boundary)
            if inter is None or inter.isEmpty():
                return None
            best = None
            for v in inter.vertices():
                pt = QgsPointXY(v.x(), v.y())
                d = origin.distance(pt)
                if 1e-6 < d <= d_cap + 1e-6 and (best is None or d < best[0]):
                    best = (d, pt)
            return best

        def sample_circle(origin: QgsPointXY, z_origin: float):
            samples = []
            deg = 0.0
            while deg < 360.0 - 1e-12:
                bearing = deg
                pt = self._point_at(origin, bearing, chord_length)
                z = sample_z(pt)
                if z is not None and z < z_origin:
                    grade = (z_origin - z) / chord_length
                    samples.append((bearing, pt, z, grade))
                deg += arc_step_deg
            return samples

        def orientation_bearing(origin: QgsPointXY, samples: list):
            """Recomputed steepest-down; ties → angular mid. Stored on node for audit."""
            if not samples:
                return None
            max_g = max(s[3] for s in samples)
            eps = 1e-12
            tied = [s for s in samples if abs(s[3] - max_g) <= eps]
            if len(tied) == 1:
                return tied[0][0]
            mid = self._angular_midpoint_deg([s[0] for s in tied])
            _ = sample_z(self._point_at(origin, mid, chord_length))
            return mid

        def chord_grade_pick(origin: QgsPointXY, z_origin: float, bearings: list[float]):
            best = None
            for bearing in bearings:
                sample = self._point_at(origin, bearing, chord_length)
                z = sample_z(sample)
                if z is None or z >= z_origin:
                    continue
                grade = (z_origin - z) / chord_length
                diff = abs(grade - g_target)
                cand = (diff, grade, bearing, sample, z)
                if best is None or cand[0] < best[0]:
                    best = cand
            return best

        def bearings_fan(orient_deg: float, span_deg: float, sidx: int):
            fan = []
            deg = 0.0
            sign = -1.0 if sidx == 0 else 1.0
            while deg < span_deg - 1e-12:
                fan.append((orient_deg + sign * deg + 360.0) % 360.0)
                deg += arc_step_deg
            fan.append((orient_deg + sign * span_deg + 360.0) % 360.0)
            return fan

        def bearings_full_circle():
            out = []
            deg = 0.0
            while deg < 360.0 - 1e-12:
                out.append(deg)
                deg += arc_step_deg
            return out

        def land_or_edge(origin: QgsPointXY, z_origin: float, pick):
            _diff, g_chord, bearing, sample, z_chord = pick
            if self._segment_inside(corridor, origin, sample):
                return sample, z_chord, g_chord, False
            hit = edge_hit(origin, bearing, chord_length)
            if hit is None:
                return None
            d, pt = hit
            z = sample_z(pt)
            if z is None or z >= z_origin:
                return None
            return pt, z, (z_origin - z) / d, True

        def edge_finale(origin: QgsPointXY, z_origin: float):
            pick = chord_grade_pick(origin, z_origin, bearings_full_circle())
            if pick is None:
                return None
            _diff, g_chord, bearing, sample, z_chord = pick
            hit = edge_hit(origin, bearing, chord_length)
            if hit is None:
                if self._segment_inside(corridor, origin, sample):
                    return sample, z_chord, g_chord, True
                return None
            d, pt = hit
            z = sample_z(pt)
            if z is None or z >= z_origin:
                return None
            return pt, z, (z_origin - z) / d, True

        def debug_rays_for_node(origin: QgsPointXY, z_origin: float, orient_deg: float, sidx: int, pick):
            """Single fall-line ray + full 180° pick fan (mark winning bearing)."""
            rays = []
            fl_pt = self._point_at(origin, orient_deg, chord_length)
            z_fl = sample_z(fl_pt)
            if z_fl is None:
                raise QgsProcessingException('DEM NA at fall-line sample point.')
            g_fl = ((z_origin - z_fl) / chord_length) * 100.0
            rays.append((
                QgsPointXY(origin),
                QgsPointXY(fl_pt),
                float(orient_deg),
                float(g_fl),
                1,  # is_fallline
                0,  # is_fan
                0,  # is_picked
            ))
            picked_bearing = None if pick is None else pick[2]
            for bearing in bearings_fan(orient_deg, 180.0, sidx):
                pt = self._point_at(origin, bearing, chord_length)
                z = sample_z(pt)
                if z is None:
                    continue
                grade_pct = ((z_origin - z) / chord_length) * 100.0
                is_picked = (
                    picked_bearing is not None
                    and abs(((bearing - picked_bearing + 180.0) % 360.0) - 180.0) < 1e-6
                )
                rays.append((
                    QgsPointXY(origin),
                    QgsPointXY(pt),
                    float(bearing),
                    float(grade_pct),
                    0,
                    1,
                    int(1 if is_picked else 0),
                ))
            return rays

        def next_node(origin: QgsPointXY, z_origin: float, sidx: int):
            """
            Returns (pt, z, grade, at_edge, orient_deg) or None.
            Orientation is always recomputed at origin.
            Stores fall-line + 180° fan rays on the current node.
            """
            samples = sample_circle(origin, z_origin)
            orient = orientation_bearing(origin, samples)
            node_i = len(nodes) - 1
            if orient is None:
                nodes[node_i]['debug_rays'] = []
                fin = edge_finale(origin, z_origin)
                if fin is None:
                    return None
                pt, z, g, edge = fin
                return pt, z, g, edge, orient

            pick = chord_grade_pick(origin, z_origin, bearings_fan(orient, 180.0, sidx))
            nodes[node_i]['debug_rays'] = debug_rays_for_node(
                origin, z_origin, orient, sidx, pick
            )
            if pick is not None:
                landed = land_or_edge(origin, z_origin, pick)
                if landed is not None:
                    pt, z, g, edge = landed
                    return pt, z, g, edge, orient
            fin = edge_finale(origin, z_origin)
            if fin is None:
                return None
            pt, z, g, edge = fin
            return pt, z, g, edge, orient

        def set_side(sidx: int):
            nonlocal side_idx, side_label
            side_idx = sidx
            side_label = 'L' if side_idx == 0 else 'R'

        def flip_side():
            set_side(1 - side_idx)

        z_start = sample_z(start_xy)
        if z_start is None:
            fail('DEM sample is NA/nodata at start point.')

        # Per node: xy, z, grade_into, at_edge, length_m, side_into,
        #           orient_deg (on departure), exit_side (0/1 on departure)
        nodes = [{
            'xy': start_xy,
            'z': z_start,
            'grade': None,
            'at_edge': False,
            'length_m': None,
            'side_into': side_label,
            'orient_deg': None,
            'exit_side': None,
            'debug_rays': [],
        }]
        exhausted = set()  # (node_index, side_idx) already tried as exit
        cur = start_xy
        z_cur = z_start
        consecutive_short = 0
        max_steps = 10000
        tip_drop = 4
        recorded_tries = []  # list of (verts, z0, z1, reaches, reason, debug_rays)

        def snapshot_debug_rays():
            out = []
            for node_i, node in enumerate(nodes):
                for ray in node.get('debug_rays') or []:
                    out.append((int(node_i),) + tuple(ray))
            return out

        def record_try(reason: str, reaches: bool = False):
            if len(nodes) < 2:
                return
            verts = [QgsPointXY(n['xy']) for n in nodes]
            recorded_tries.append((
                verts,
                float(nodes[0]['z']),
                float(nodes[-1]['z']),
                reaches,
                reason,
                snapshot_debug_rays(),
            ))

        def reaches_end(xy: QgsPointXY) -> bool:
            return xy.distance(end_xy) <= end_buffer + 1e-9

        def try_backtrack(reason: str) -> bool:
            """
            Record current try, drop failed tip, rewind to an interior node,
            try opposite exit side once per (node, side).
            Returns False when rewind exhausts back at start (terminal stop).
            """
            nonlocal cur, z_cur, consecutive_short
            feedback.pushInfo(f'Backtrack ({reason}).')
            record_try(reason, reaches=False)

            if len(nodes) <= tip_drop:
                feedback.pushInfo(
                    'End: rewind blocked (path too short to drop tip) — at/near start.'
                )
                return False

            del nodes[-tip_drop:]

            while len(nodes) >= 1:
                i = len(nodes) - 1
                node = nodes[i]

                if node['at_edge'] and len(nodes) > 1:
                    nodes.pop()
                    continue

                prev_exit = node.get('exit_side')
                candidates = []
                if prev_exit is not None:
                    candidates.append(1 - prev_exit)
                    candidates.append(prev_exit)
                else:
                    candidates.append(1 - side_idx)
                    candidates.append(side_idx)

                chosen = None
                for s in candidates:
                    if (i, s) not in exhausted:
                        chosen = s
                        break

                if chosen is None:
                    if i == 0 or len(nodes) == 1:
                        feedback.pushInfo(
                            'End: rewind returned to start — all opening sides exhausted.'
                        )
                        return False
                    nodes.pop()
                    continue

                set_side(chosen)
                cur = node['xy']
                z_cur = node['z']
                node['exit_side'] = None
                consecutive_short = 0
                feedback.pushInfo(
                    f'Rewind to node {i}, try side {side_label} '
                    f'(stored orient={node.get("orient_deg")}).'
                )
                return True

            feedback.pushInfo('End: rewind emptied path back to start.')
            return False

        reached_end = False

        for step in range(max_steps):
            if feedback.isCanceled():
                break

            if reaches_end(cur):
                reached_end = True
                feedback.pushInfo(
                    f'End reached within buffer {end_buffer:.1f} m '
                    f'(dist={cur.distance(end_xy):.2f} m).'
                )
                record_try('reached_end', reaches=True)
                break

            parent_i = len(nodes) - 1
            nxt = next_node(cur, z_cur, side_idx)
            if nxt is None:
                if try_backtrack('no candidate'):
                    continue
                record_try('no_candidate_final', reaches=False)
                break

            pt, z_pt, grade, at_edge, orient = nxt
            length_m = cur.distance(pt)

            nodes[parent_i]['orient_deg'] = orient
            nodes[parent_i]['exit_side'] = side_idx
            exhausted.add((parent_i, side_idx))

            nodes.append({
                'xy': pt,
                'z': z_pt,
                'grade': grade,
                'at_edge': at_edge,
                'length_m': length_m,
                'side_into': side_label,
                'orient_deg': None,
                'exit_side': None,
                'debug_rays': [],
            })

            cur = pt
            z_cur = z_pt

            if reaches_end(cur):
                reached_end = True
                feedback.pushInfo(
                    f'End reached within buffer {end_buffer:.1f} m '
                    f'(dist={cur.distance(end_xy):.2f} m).'
                )
                record_try('reached_end', reaches=True)
                break

            if length_m < chord_length - 1e-9:
                consecutive_short += 1
                if consecutive_short >= 2:
                    if try_backtrack('two short chords'):
                        continue
                    record_try('two_short_final', reaches=False)
                    break
            else:
                consecutive_short = 0

            if at_edge:
                flip_side()
                feedback.pushInfo(
                    f'Edge node {len(nodes) - 1}: flipped to side {side_label}, continuing.'
                )
        else:
            feedback.pushWarning('Hit max steps.')
            record_try('max_steps', reaches=False)

        # Final path if not already recorded as a try
        if recorded_tries:
            last_verts = recorded_tries[-1][0]
            cur_verts = [n['xy'] for n in nodes]
            same = (
                len(last_verts) == len(cur_verts)
                and all(a.distance(b) < 1e-6 for a, b in zip(last_verts, cur_verts))
            )
            if not same and len(nodes) >= 2:
                record_try('final', reaches=reached_end)
        elif len(nodes) >= 2:
            record_try('final', reaches=reached_end)

        results = {}

        fields = QgsFields()
        fields.append(QgsField('node_id', QVariant.Int))
        fields.append(QgsField('side', QVariant.String))
        fields.append(QgsField('exit_side', QVariant.String))
        fields.append(QgsField('length_m', QVariant.Double))
        fields.append(QgsField('orient_deg', QVariant.Double))
        fields.append(QgsField('grade_pct', QVariant.Double))
        fields.append(QgsField('target_pct', QVariant.Double))
        fields.append(QgsField('grade_err', QVariant.Double))
        fields.append(QgsField('at_edge', QVariant.Int))
        fields.append(QgsField('reaches_end', QVariant.Int))

        sink, sink_id = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            fields,
            QgsWkbTypes.Point,
            corridor_crs,
        )
        if sink is not None:
            for i, node in enumerate(nodes):
                grade = node['grade']
                if grade is None:
                    g_pct = None
                    g_err = None
                else:
                    g_pct = grade * 100.0
                    g_err = g_pct - grade_pct
                exit_side = node.get('exit_side')
                exit_side_str = None if exit_side is None else ('L' if exit_side == 0 else 'R')
                feat = QgsFeature(fields)
                feat.setGeometry(QgsGeometry.fromPointXY(node['xy']))
                feat.setAttributes([
                    int(i),
                    node['side_into'],
                    exit_side_str,
                    node['length_m'],
                    node.get('orient_deg'),
                    g_pct,
                    float(grade_pct),
                    g_err,
                    int(1 if node['at_edge'] else 0),
                    int(1 if reached_end else 0),
                ])
                sink.addFeature(feat, QgsFeatureSink.FastInsert)
            results[self.OUTPUT] = sink_id

        best_try_idx = None
        best_try_len = -1.0
        best_end_idx = None
        try_path_lens = []
        for try_id, (verts, z0, z1, reaches, reason, _rays) in enumerate(recorded_tries):
            if len(verts) < 2:
                try_path_lens.append(None)
                continue
            path_len = sum(verts[i].distance(verts[i + 1]) for i in range(len(verts) - 1))
            try_path_lens.append((path_len, z0, z1, reaches, reason, verts))
            if reaches:
                best_end_idx = try_id
            if path_len > best_try_len:
                best_try_len = path_len
                best_try_idx = try_id

        line_fields = QgsFields()
        line_fields.append(QgsField('try_id', QVariant.Int))
        line_fields.append(QgsField('n_seg', QVariant.Int))
        line_fields.append(QgsField('length_m', QVariant.Double))
        line_fields.append(QgsField('avg_grade_pct', QVariant.Double))
        line_fields.append(QgsField('reaches_end', QVariant.Int))
        line_fields.append(QgsField('reason', QVariant.String))

        line_sink, line_id = self.parameterAsSink(
            parameters,
            self.LINE_OUTPUT,
            context,
            line_fields,
            QgsWkbTypes.LineString,
            corridor_crs,
        )
        if line_sink is not None:
            n_tries_out = 0
            for try_id, packed in enumerate(try_path_lens):
                if packed is None:
                    continue
                path_len, z0, z1, reaches, reason, verts = packed
                avg_grade_pct = ((z0 - z1) / path_len) * 100.0
                line_feat = QgsFeature(line_fields)
                line_feat.setGeometry(QgsGeometry.fromPolylineXY(verts))
                line_feat.setAttributes([
                    int(try_id),
                    int(len(verts) - 1),
                    float(path_len),
                    avg_grade_pct,
                    int(1 if reaches else 0),
                    reason,
                ])
                line_sink.addFeature(line_feat, QgsFeatureSink.FastInsert)
                n_tries_out += 1
            feedback.pushInfo(f'Path tries written: {n_tries_out}')
            results[self.LINE_OUTPUT] = line_id

            _style_path_tries_layer(_get_output_layer(context, line_id))
            if context.project():
                details = QgsProcessingContext.LayerDetails(
                    self.tr('Path tries'),
                    context.project(),
                    self.LINE_OUTPUT,
                )
                details.setPostProcessor(_PathTriesStylePostProcessor())
                context.addLayerToLoadOnCompletion(line_id, details)

        fl_fields = QgsFields()
        fl_fields.append(QgsField('node_id', QVariant.Int))
        fl_fields.append(QgsField('bearing', QVariant.Double))
        fl_fields.append(QgsField('grade_pct', QVariant.Double))
        fl_fields.append(QgsField('is_fallline', QVariant.Int))
        fl_fields.append(QgsField('is_fan', QVariant.Int))
        fl_fields.append(QgsField('is_picked', QVariant.Int))
        fl_fields.append(QgsField('try_id', QVariant.Int))

        fl_sink, fl_id = self.parameterAsSink(
            parameters,
            self.FALLLINE_OUTPUT,
            context,
            fl_fields,
            QgsWkbTypes.LineString,
            corridor_crs,
        )
        if fl_sink is not None:
            debug_try_id = best_end_idx if best_end_idx is not None else best_try_idx
            n_fl = 0
            if debug_try_id is not None:
                _v, _z0, _z1, _r, _reason, rays = recorded_tries[debug_try_id]
                for node_i, a, b, bearing, grade_pct, is_fl, is_fan, is_picked in rays:
                    feat = QgsFeature(fl_fields)
                    feat.setGeometry(QgsGeometry.fromPolylineXY([a, b]))
                    feat.setAttributes([
                        int(node_i),
                        float(bearing),
                        float(grade_pct),
                        int(is_fl),
                        int(is_fan),
                        int(is_picked),
                        int(debug_try_id),
                    ])
                    fl_sink.addFeature(feat, QgsFeatureSink.FastInsert)
                    n_fl += 1
                which = 'end-reached' if best_end_idx is not None else 'longest'
                feedback.pushInfo(
                    f'Debug rays written: {n_fl} ({which} try_id={debug_try_id})'
                )
            else:
                feedback.pushInfo('Debug rays written: 0 (no try path)')
            results[self.FALLLINE_OUTPUT] = fl_id

        if reached_end:
            feedback.pushInfo(f'Wrote {len(nodes)} node(s); end reached.')
        else:
            feedback.pushWarning(
                f'Wrote {len(nodes)} node(s); end not reached '
                f'(dist tip→end={nodes[-1]["xy"].distance(end_xy):.2f} m).'
            )
        return results
