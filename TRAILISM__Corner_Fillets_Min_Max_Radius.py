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

from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtGui import QColor
from qgis.core import (
    QgsProcessing,
    QgsFeatureSink,
    QgsProcessingException,
    QgsProcessingAlgorithm,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFeatureSink,
    QgsGeometry,
    QgsPoint,
    QgsPointXY,
    QgsLineString,
    QgsGeometryUtils,
    QgsCompoundCurve,
    QgsCircularString,
    QgsTriangle,
    QgsFeature,
    QgsLineSymbol,
    QgsSimpleLineSymbolLayer,
    QgsSingleSymbolRenderer,
    QgsUnitTypes,
)
from math import sin, tan, pi, sqrt, atan2, asin, ceil, cos


class CornerFilletsMinMaxRadiusAlgorithm(QgsProcessingAlgorithm):
    """
    Classical **corner fillets** (not vertex-centered circles):
    - Optional cleaning: minimum spacing between input vertices along the line.
    - At each corner, the fillet center lies on the **interior angle bisector**
      (via triangle **incenter** direction from the vertex).
    - The circular arc is **tangent to the two incident straight legs** (tangent
      points lie on those segments).
    - Radius is derived from local geometry and clamped by user min/max radius;
      arcs are linearized to a max chord length.

    For **vertex-centered circles** with straight **outer common tangents**, use
    the separate algorithm TRAILISM: Corner Arcs (vertex circles, outside).
    """

    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'
    MIN_RADIUS = 'MIN_RADIUS'
    MAX_RADIUS = 'MAX_RADIUS'
    MAX_SEG_LEN = 'MAX_SEG_LEN'
    MIN_VERTEX_SPACING = 'MIN_VERTEX_SPACING'
    DRAW_BERMS = 'DRAW_BERMS'
    OUTPUT_BERMS = 'OUTPUT_BERMS'
    BERM_OFFSET = 'BERM_OFFSET'
    BERM_RADIUS_THRESHOLD = 'BERM_RADIUS_THRESHOLD'
    MIN_BERM_LENGTH = 'MIN_BERM_LENGTH'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return CornerFilletsMinMaxRadiusAlgorithm()

    def name(self):
        return 'corner_fillets_min_max_radius'

    def displayName(self):
        return self.tr('TRAILISM: Corner Fillets (min/max radius)')

    def shortHelpString(self):
        return self.tr(
            'TRAILISM ALGORITHM (bisector fillets):\n\n'
            '- Optional: cleans vertices that are closer than a minimum spacing along the polyline.\n'
            '- At each corner, builds a **fillet**: the arc lies inside the turn; its center sits on the '
            '**interior angle bisector** (incenter-based), and the arc is **tangent to both legs** '
            '(tangent points on the adjacent segments).\n'
            '- Fillet radius uses the same geometric allowance as “Round corner by min/max radius”, '
            'then clamps to your Min/Max radius.\n'
            '- Arcs are approximated with straight segments so chord length ≤ Max chord length.\n'
            '- Optional berm output: offset arc outside the fillet (same as the original round-corner tool).\n\n'
            'This is **not** vertex-centered rounding. For circles centered on each vertex with outer tangents, '
            'use “TRAILISM: Corner Arcs (vertex circles, outside)”.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorLine]
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MIN_VERTEX_SPACING,
                self.tr('Min vertex spacing before filleting (map units)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MIN_RADIUS,
                self.tr('Min radius (0 = no minimum)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=3.0,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_RADIUS,
                self.tr('Max radius (0 = no maximum)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=6.0,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_SEG_LEN,
                self.tr('Max chord length for arcs (0 = default tolerance)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.DRAW_BERMS,
                self.tr('Draw berms as offset arc (for curves with radius ≤ threshold)'),
                defaultValue=False
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.BERM_OFFSET,
                self.tr('Berm offset distance (map units)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.BERM_RADIUS_THRESHOLD,
                self.tr('Berm radius threshold (draw berm when fillet radius ≤ value)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=6.0,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MIN_BERM_LENGTH,
                self.tr('Min berm length (map units); only output berms ≥ this length'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=5.0,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_BERMS,
                self.tr('Berms layer (optional)')
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Curved lines (filleted)')
            )
        )

    def proper_angle(self, angle: float) -> float:
        """
        Constrain angle to [0, pi] radians (use interior angle).
        """
        return angle if 0 <= angle <= pi else 2 * pi - angle

    def _segment_length(self, a: QgsPointXY, b: QgsPointXY) -> float:
        """
        Euclidean length of segment AB in 2D.
        """
        return sqrt(QgsGeometryUtils.sqrDistance2D(QgsPoint(a), QgsPoint(b)))

    def _segment_param(self, a: QgsPointXY, b: QgsPointXY, p) -> float:
        """
        Scalar projection parameter s of point P onto segment AB where
        P ≈ A + s*(B-A). Returns s clamped to [0, 1]. Accepts P as QgsPoint or QgsPointXY.
        """
        ax, ay = a.x(), a.y()
        bx, by = b.x(), b.y()
        px, py = (p.x(), p.y())
        vx, vy = bx - ax, by - ay
        wx, wy = px - ax, py - ay
        denom = vx * vx + vy * vy
        if denom <= 0.0:
            return 0.0
        s = (wx * vx + wy * vy) / denom
        if s < 0.0:
            return 0.0
        if s > 1.0:
            return 1.0
        return s

    def clean_vertices_min_spacing(self, points: list, min_spacing: float) -> list:
        """
        Simple cleaning pass: remove vertices that are closer than min_spacing
        to the last kept vertex along the polyline.
        """
        if len(points) <= 2 or min_spacing <= 0.0:
            return points

        cleaned = [points[0]]
        last_keep = points[0]
        for p in points[1:]:
            if self._segment_length(last_keep, p) >= min_spacing:
                cleaned.append(p)
                last_keep = p
        if cleaned[-1] != points[-1]:
            cleaned.append(points[-1])
        return cleaned

    def preprocess_remove_problematic_vertices(self, points: list, min_radius: float, max_radius: float) -> list:
        """
        Same as in the max-radius algorithm: remove vertices where arcs cannot fit
        because of blunt angles combined with short adjacent segments.
        """
        if len(points) < 3:
            return points

        filtered = [points[0]]  # Always keep first point
        i = 1

        while i < len(points) - 1:
            p_prev = points[i - 1]
            p_mid = points[i]
            p_next = points[i + 1]

            theta = self.proper_angle(
                QgsGeometryUtils.angleBetweenThreePoints(
                    p_prev.x(), p_prev.y(),
                    p_mid.x(), p_mid.y(),
                    p_next.x(), p_next.y()
                )
            )

            len_prev = self._segment_length(p_prev, p_mid)
            len_next = self._segment_length(p_mid, p_next)

            if len_prev < 1e-9 or len_next < 1e-9:
                i += 1
                continue

            tan_half = tan(theta / 2.0)
            if tan_half <= 1e-12:
                i += 1
                continue

            r_allowed = 0.5 * min(len_prev, len_next) * tan_half
            if r_allowed <= 0.0:
                i += 1
                continue

            r_cap = r_allowed
            if max_radius and max_radius > 0.0:
                r_cap = min(r_cap, max_radius)
            if min_radius and min_radius > 0.0 and r_cap < min_radius:
                r_estimate = min_radius
            else:
                r_estimate = r_cap

            t_prev_required = r_estimate / tan_half
            t_next_required = r_estimate / tan_half

            if t_prev_required >= len_prev or t_next_required >= len_next:
                i += 1
                continue

            filtered.append(p_mid)
            i += 1

        filtered.append(points[-1])
        return filtered

    def _segment_arc_by_chord(self, center_pt: QgsPoint, t1: QgsPoint, t2: QgsPoint, phi_expected: float, max_chord_len: float) -> QgsLineString:
        """
        Segment a circular arc defined by center and two points t1->t2
        into a linestring with max chord length <= max_chord_len.
        """
        cx, cy = center_pt.x(), center_pt.y()
        a1 = atan2(t1.y() - cy, t1.x() - cx)
        a2 = atan2(t2.y() - cy, t2.x() - cx)

        delta = a2 - a1
        while delta <= -pi:
            delta += 2 * pi
        while delta > pi:
            delta -= 2 * pi

        if abs(abs(delta) - phi_expected) > abs(abs(-delta) - phi_expected):
            delta = -delta

        R = sqrt((t1.x() - cx) ** 2 + (t1.y() - cy) ** 2)
        if R <= 0.0 or max_chord_len <= 0.0:
            return QgsLineString([t1, t2])

        step = 2.0 * asin(min(1.0, max_chord_len / (2.0 * max(R, 1e-12))))
        if step <= 0.0:
            return QgsLineString([t1, t2])

        n = max(1, int(ceil(abs(delta) / step)))
        step_signed = delta / n

        pts = [t1]
        for k in range(1, n):
            ang = a1 + step_signed * k
            x = cx + R * cos(ang)
            y = cy + R * sin(ang)
            pts.append(QgsPoint(x, y))
        pts.append(t2)

        return QgsLineString(pts)

    def smooth_line_corners(self, points: list, max_seg_len: float, min_radius: float, max_radius: float, draw_berms: bool, berm_geoms: list, berm_offset: float = 1.0, berm_radius_threshold: float = 6.0) -> QgsGeometry:
        """
        Bisector-based fillet construction (see RoundCornerMaxRadiusAlgorithm),
        after optional vertex spacing clean and preprocess_remove_problematic_vertices.
        """
        if len(points) < 3:
            return QgsGeometry.fromPolyline(list(map(lambda p: QgsPoint(p), points)))

        points = self.preprocess_remove_problematic_vertices(points, min_radius, max_radius)

        if len(points) < 3:
            return QgsGeometry.fromPolyline(list(map(lambda p: QgsPoint(p), points)))

        compound = QgsCompoundCurve()
        current_point = QgsPoint(points[0])
        compound.addVertex(current_point)

        for i in range(1, len(points) - 1):
            p_prev = points[i - 1]
            p_mid = points[i]
            p_next = points[i + 1]

            theta = self.proper_angle(
                QgsGeometryUtils.angleBetweenThreePoints(
                    p_prev.x(), p_prev.y(),
                    p_mid.x(), p_mid.y(),
                    p_next.x(), p_next.y()
                )
            )

            len_prev = self._segment_length(p_prev, p_mid)
            len_next = self._segment_length(p_mid, p_next)

            if len_prev < 1e-9 or len_next < 1e-9 or abs(tan(theta / 2.0)) < 1e-12:
                line = QgsLineString([current_point, QgsPoint(p_mid)])
                compound.addCurve(line, extendPrevious=True)
                current_point = QgsPoint(p_mid)
                continue

            tan_half = tan(theta / 2.0)
            if tan_half <= 1e-12:
                line = QgsLineString([current_point, QgsPoint(p_mid)])
                compound.addCurve(line, extendPrevious=True)
                current_point = QgsPoint(p_mid)
                continue
            r_allowed = 0.5 * min(len_prev, len_next) * tan_half

            if r_allowed <= 0.0:
                line = QgsLineString([current_point, QgsPoint(p_mid)])
                compound.addCurve(line, extendPrevious=True)
                current_point = QgsPoint(p_mid)
                continue

            r_cap = r_allowed
            if max_radius and max_radius > 0.0:
                r_cap = min(r_cap, max_radius)
            if min_radius and min_radius > 0.0 and r_cap < min_radius:
                r_final = min_radius
            else:
                r_final = r_cap

            triangle = QgsTriangle(p_prev, p_mid, p_next)
            incenter = QgsPointXY(triangle.inscribedCircle().center())
            d_center = r_final / max(sin(theta / 2.0), 1e-12)

            center_A = QgsGeometryUtils.pointOnLineWithDistance(QgsPoint(p_mid), QgsPoint(incenter), d_center)
            t1_A = QgsGeometryUtils.projectPointOnSegment(center_A, QgsPoint(p_prev), QgsPoint(p_mid))
            t2_A = QgsGeometryUtils.projectPointOnSegment(center_A, QgsPoint(p_mid), QgsPoint(p_next))

            center_B = QgsGeometryUtils.pointOnLineWithDistance(QgsPoint(p_mid), QgsPoint(incenter), -d_center)
            t1_B = QgsGeometryUtils.projectPointOnSegment(center_B, QgsPoint(p_prev), QgsPoint(p_mid))
            t2_B = QgsGeometryUtils.projectPointOnSegment(center_B, QgsPoint(p_mid), QgsPoint(p_next))

            def central_angle(center_point, a_point, b_point):
                ang = QgsGeometryUtils.angleBetweenThreePoints(
                    a_point.x(), a_point.y(),
                    center_point.x(), center_point.y(),
                    b_point.x(), b_point.y()
                )
                return self.proper_angle(ang)

            phi_expected = pi - theta
            phi_A = central_angle(center_A, t1_A, t2_A)
            phi_B = central_angle(center_B, t1_B, t2_B)

            if abs(phi_A - phi_expected) <= abs(phi_B - phi_expected):
                center_pt, t1, t2 = center_A, t1_A, t2_A
            else:
                center_pt, t1, t2 = center_B, t1_B, t2_B

            if min_radius and min_radius > 0.0 and r_final >= min_radius and (r_allowed < min_radius):
                def clamp_to_segment(a_xy: QgsPointXY, b_xy: QgsPointXY, p: QgsPoint) -> QgsPoint:
                    s = self._segment_param(a_xy, b_xy, p)
                    if s <= 0.0:
                        return QgsPoint(a_xy)
                    if s >= 1.0:
                        return QgsPoint(b_xy)
                    ax, ay = a_xy.x(), a_xy.y()
                    bx, by = b_xy.x(), b_xy.y()
                    return QgsPoint(ax + (bx - ax) * s, ay + (by - ay) * s)

                t1 = clamp_to_segment(p_prev, p_mid, t1)
                t2 = clamp_to_segment(p_mid, p_next, t2)

            if QgsGeometryUtils.sqrDistance2D(current_point, t1) > 1e-18:
                line_to_t1 = QgsLineString([current_point, t1])
                compound.addCurve(line_to_t1, extendPrevious=True)

            if max_seg_len and max_seg_len > 0.0:
                arc_line = self._segment_arc_by_chord(center_pt, t1, t2, phi_expected, max_seg_len)
                compound.addCurve(arc_line, extendPrevious=True)
            else:
                arc = QgsCircularString.fromTwoPointsAndCenter(t1, t2, center_pt, useShortestArc=True)
                compound.addCurve(arc.curveToLine(0.1), extendPrevious=True)

            if draw_berms and r_final <= berm_radius_threshold:
                cx, cy = center_pt.x(), center_pt.y()
                a1 = atan2(t1.y() - cy, t1.x() - cx)
                a2 = atan2(t2.y() - cy, t2.x() - cx)
                R_off = r_final + berm_offset
                t1_off = QgsPoint(cx + R_off * cos(a1), cy + R_off * sin(a1))
                t2_off = QgsPoint(cx + R_off * cos(a2), cy + R_off * sin(a2))

                if max_seg_len and max_seg_len > 0.0:
                    berm_line = self._segment_arc_by_chord(center_pt, t1_off, t2_off, phi_expected, max_seg_len)
                else:
                    berm_arc = QgsCircularString.fromTwoPointsAndCenter(t1_off, t2_off, center_pt, useShortestArc=True)
                    berm_line = berm_arc.curveToLine(0.1)
                berm_geoms.append(berm_line)

            current_point = t2

        if QgsGeometryUtils.sqrDistance2D(current_point, QgsPoint(points[-1])) > 1e-18:
            end_line = QgsLineString([current_point, QgsPoint(points[-1])])
            compound.addCurve(end_line, extendPrevious=True)

        return compound.curveToLine()

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        max_seg_len = self.parameterAsDouble(parameters, self.MAX_SEG_LEN, context)
        min_radius = self.parameterAsDouble(parameters, self.MIN_RADIUS, context)
        max_radius = self.parameterAsDouble(parameters, self.MAX_RADIUS, context)
        min_vertex_spacing = self.parameterAsDouble(parameters, self.MIN_VERTEX_SPACING, context)

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            source.fields(),
            source.wkbType(),
            source.sourceCrs()
        )

        draw_berms = self.parameterAsBool(parameters, self.DRAW_BERMS, context)
        berm_offset = self.parameterAsDouble(parameters, self.BERM_OFFSET, context)
        berm_radius_threshold = self.parameterAsDouble(parameters, self.BERM_RADIUS_THRESHOLD, context)
        min_berm_length = self.parameterAsDouble(parameters, self.MIN_BERM_LENGTH, context)
        min_berm_length = float(min_berm_length) if min_berm_length is not None else 5.0
        berm_sink = None
        berm_dest_id = None
        if draw_berms:
            (berm_sink, berm_dest_id) = self.parameterAsSink(
                parameters,
                self.OUTPUT_BERMS,
                context,
                source.fields(),
                source.wkbType(),
                source.sourceCrs()
            )

        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))

        total = 100.0 / source.featureCount() if source.featureCount() else 0

        for current, feature in enumerate(source.getFeatures()):
            if feedback.isCanceled():
                break

            geom = feature.geometry()
            if not geom:
                continue

            pts = [QgsPointXY(v.x(), v.y()) for v in geom.vertices()]

            if len(pts) < 3:
                new_geom = QgsGeometry.fromPolyline(list(map(lambda p: QgsPoint(p), pts)))
            else:
                cleaned = self.clean_vertices_min_spacing(pts, min_vertex_spacing)
                if len(cleaned) < 3:
                    new_geom = QgsGeometry.fromPolyline(list(map(lambda p: QgsPoint(p), cleaned)))
                    berm_geoms = []
                else:
                    berm_geoms = []
                    new_geom = self.smooth_line_corners(cleaned, max_seg_len, min_radius, max_radius, draw_berms, berm_geoms, berm_offset=berm_offset, berm_radius_threshold=berm_radius_threshold)

            new_feature = QgsFeature(feature.fields())
            new_feature.setAttributes(feature.attributes())
            new_feature.setGeometry(new_geom)
            sink.addFeature(new_feature, QgsFeatureSink.FastInsert)

            if draw_berms and berm_sink is not None and berm_geoms:
                for bl in berm_geoms:
                    berm_geom = QgsGeometry(bl)
                    if berm_geom.length() < min_berm_length:
                        continue
                    berm_feat = QgsFeature(feature.fields())
                    berm_feat.setAttributes(feature.attributes())
                    berm_feat.setGeometry(berm_geom)
                    berm_sink.addFeature(berm_feat, QgsFeatureSink.FastInsert)

            feedback.setProgress(int(current * total))

        layer = context.getMapLayer(dest_id)
        if layer is not None:
            symbol_layer = QgsSimpleLineSymbolLayer()
            symbol_layer.setColor(QColor(255, 255, 255))
            symbol_layer.setWidth(0.3)
            symbol_layer.setWidthUnit(QgsUnitTypes.RenderMillimeters)

            symbol = QgsLineSymbol()
            symbol.changeSymbolLayer(0, symbol_layer)
            renderer = QgsSingleSymbolRenderer(symbol)
            layer.setRenderer(renderer)
            layer.triggerRepaint()

        results = {self.OUTPUT: dest_id}
        if draw_berms and berm_dest_id is not None:
            results[self.OUTPUT_BERMS] = berm_dest_id
            berm_layer = context.getMapLayer(berm_dest_id)
            if berm_layer is not None:
                b_sym_layer = QgsSimpleLineSymbolLayer()
                b_sym_layer.setColor(QColor(255, 255, 255))
                b_sym_layer.setWidth(0.5)
                b_sym_layer.setWidthUnit(QgsUnitTypes.RenderMillimeters)

                b_symbol = QgsLineSymbol()
                b_symbol.changeSymbolLayer(0, b_sym_layer)
                b_renderer = QgsSingleSymbolRenderer(b_symbol)
                berm_layer.setRenderer(b_renderer)
                berm_layer.triggerRepaint()
        return results
