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


class RoundCornerFixedRadiusAlgorithm(QgsProcessingAlgorithm):
    """
    Inserts circular arcs at each vertex of a polyline using a FIXED radius.
    The radius is applied to all corners, trimmed to fit within adjacent segments.
    """

    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'
    RADIUS = 'RADIUS'
    MAX_SEG_LEN = 'MAX_SEG_LEN'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return RoundCornerFixedRadiusAlgorithm()

    def name(self):
        return 'roundcornerfixedradius'

    def displayName(self):
        return self.tr('TRAILISM: Round Corner with Fixed Radius')

    def shortHelpString(self):
        return self.tr(
            'TRAILISM ALGORITHM:\n\n'
            'Rounds corners with a fixed radius applied to all vertices.\n\n'
            'Key terms:\n'
            '- Fixed radius R: fillet radius applied to every corner.\n'
            '- Sweep angle φ: inside arc angle, φ = π − θ, where θ is the interior angle at the vertex.\n'
            '- Chord length c: straight segment between two neighboring points on the arc.\n'
            '- Deviation tolerance t: max distance between true arc and line approximation.\n\n'
            'Parameters:\n'
            '- Fixed radius: radius applied to all corners (map units).\n'
            '  If the radius is too large for a corner, it will be trimmed to fit within the adjacent segments.\n'
            '- Max chord length for arcs (0 = default tolerance): if > 0, arc is segmented so each chord ≤ given length (map units).\n'
            '  If 0, a deviation-based approximation with t = 0.1 map units is used.'
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
                self.RADIUS,
                self.tr('Fixed radius'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=5.0,
                minValue=0.1
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_SEG_LEN,
                self.tr('Max chord length for arcs (0 = default tolerance)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.2,
                minValue=0.0
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Curved Lines with Fixed Radius')
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

    def _segment_arc_by_chord(self, center_pt: QgsPoint, t1: QgsPoint, t2: QgsPoint, phi_expected: float, max_chord_len: float) -> QgsLineString:
        """
        Segment a circular arc defined by center and two points t1->t2
        into a linestring with max chord length <= max_chord_len.
        Chooses the shorter sweep matching the inside angle.
        """
        cx, cy = center_pt.x(), center_pt.y()
        a1 = atan2(t1.y() - cy, t1.x() - cx)
        a2 = atan2(t2.y() - cy, t2.x() - cx)

        # Normalize to (-pi, pi]
        delta = a2 - a1
        while delta <= -pi:
            delta += 2 * pi
        while delta > pi:
            delta -= 2 * pi

        # Ensure the sweep corresponds to the intended inside arc
        if abs(abs(delta) - phi_expected) > abs(abs(-delta) - phi_expected):
            delta = -delta

        R = sqrt((t1.x() - cx) ** 2 + (t1.y() - cy) ** 2)
        if R <= 0.0 or max_chord_len <= 0.0:
            return QgsLineString([t1, t2])

        # chord c = 2 R sin(dtheta/2) <= L  => dtheta <= 2 asin(L/(2R))
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

    def smooth_line_corners_fixed(self, points: list, max_seg_len: float, fixed_radius: float) -> QgsGeometry:
        """
        Build a smoothed polyline by replacing each corner with a circular arc
        using the fixed radius, trimmed to fit within adjacent segments.
        """
        if len(points) < 3:
            return QgsGeometry.fromPolyline(list(map(lambda p: QgsPoint(p), points)))

        compound = QgsCompoundCurve()
        current_point = QgsPoint(points[0])
        compound.addVertex(current_point)

        # Iterate over internal vertices
        for i in range(1, len(points) - 1):
            p_prev = points[i - 1]
            p_mid = points[i]
            p_next = points[i + 1]

            # Angle at the corner
            theta = self.proper_angle(
                QgsGeometryUtils.angleBetweenThreePoints(
                    p_prev.x(), p_prev.y(),
                    p_mid.x(), p_mid.y(),
                    p_next.x(), p_next.y()
                )
            )

            # Segment lengths
            len_prev = self._segment_length(p_prev, p_mid)
            len_next = self._segment_length(p_mid, p_next)

            # Skip if degenerate or near-straight
            if len_prev < 1e-9 or len_next < 1e-9 or abs(tan(theta / 2.0)) < 1e-12:
                # Just draw up to the corner
                line = QgsLineString([current_point, QgsPoint(p_mid)])
                compound.addCurve(line, extendPrevious=True)
                current_point = QgsPoint(p_mid)
                continue

            # Calculate maximum allowed radius for this corner
            tan_half = tan(theta / 2.0)
            if tan_half <= 1e-12:
                line = QgsLineString([current_point, QgsPoint(p_mid)])
                compound.addCurve(line, extendPrevious=True)
                current_point = QgsPoint(p_mid)
                continue
            
            r_max_allowed = 0.5 * min(len_prev, len_next) * tan_half
            
            # Use the smaller of fixed radius or maximum allowed
            r_final = min(fixed_radius, r_max_allowed)
            
            if r_final <= 0.0:
                line = QgsLineString([current_point, QgsPoint(p_mid)])
                compound.addCurve(line, extendPrevious=True)
                current_point = QgsPoint(p_mid)
                continue

            # Find circle center along angle bisector
            triangle = QgsTriangle(p_prev, p_mid, p_next)
            incenter = QgsPointXY(triangle.inscribedCircle().center())
            d_center = r_final / max(sin(theta / 2.0), 1e-12)

            # Candidate A: towards incenter
            center_A = QgsGeometryUtils.pointOnLineWithDistance(QgsPoint(p_mid), QgsPoint(incenter), d_center)
            t1_A = QgsGeometryUtils.projectPointOnSegment(center_A, QgsPoint(p_prev), QgsPoint(p_mid))
            t2_A = QgsGeometryUtils.projectPointOnSegment(center_A, QgsPoint(p_mid), QgsPoint(p_next))

            # Candidate B: opposite direction from incenter
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

            # Connect current_point → t1 (straight), then arc t1 → t2 with given center
            if QgsGeometryUtils.sqrDistance2D(current_point, t1) > 1e-18:
                line_to_t1 = QgsLineString([current_point, t1])
                compound.addCurve(line_to_t1, extendPrevious=True)

            if max_seg_len and max_seg_len > 0.0:
                arc_line = self._segment_arc_by_chord(center_pt, t1, t2, phi_expected, max_seg_len)
                compound.addCurve(arc_line, extendPrevious=True)
            else:
                arc = QgsCircularString.fromTwoPointsAndCenter(t1, t2, center_pt, useShortestArc=True)
                # Fallback linearization using deviation tolerance t = 0.1 map units
                compound.addCurve(arc.curveToLine(0.1), extendPrevious=True)

            current_point = t2

        # Finish with straight segment to the last original point
        if QgsGeometryUtils.sqrDistance2D(current_point, QgsPoint(points[-1])) > 1e-18:
            end_line = QgsLineString([current_point, QgsPoint(points[-1])])
            compound.addCurve(end_line, extendPrevious=True)

        return compound.curveToLine()

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        fixed_radius = self.parameterAsDouble(parameters, self.RADIUS, context)
        max_seg_len = self.parameterAsDouble(parameters, self.MAX_SEG_LEN, context)

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
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

            # Collect vertices as QgsPointXY for processing
            pts = [QgsPointXY(v.x(), v.y()) for v in geom.vertices()]

            if len(pts) < 3:
                # Nothing to round
                new_geom = QgsGeometry.fromPolyline(list(map(lambda p: QgsPoint(p), pts)))
            else:
                new_geom = self.smooth_line_corners_fixed(pts, max_seg_len, fixed_radius)

            new_feature = QgsFeature(feature.fields())
            new_feature.setAttributes(feature.attributes())
            new_feature.setGeometry(new_geom)
            sink.addFeature(new_feature, QgsFeatureSink.FastInsert)

            feedback.setProgress(int(current * total))

        # Apply default white styling with 0.3 mm line width to the output layer
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

        return {self.OUTPUT: dest_id}