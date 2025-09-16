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

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsProcessing,
                       QgsWkbTypes,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterNumber,
                       QgsGeometry,
                       QgsPoint,
                       QgsPointXY,
                       QgsLineString,
                       QgsGeometryUtils,
                       QgsCompoundCurve,
                       QgsCircularString,
                       QgsTriangle,
                       QgsField,
                       QgsFeature)
from qgis import processing
from math import sin, pi


class FindCurvesAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm finds curves in a vector line layer.
    """
    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return FindCurvesAlgorithm()

    def name(self):
        return 'find curves'

    def displayName(self):
        return self.tr('TRAILISM: Find Curves')

    def shortHelpString(self):
        return self.tr(
            'TRAILISM ALGORITHM:\n\n'
            'Finds curves in a vector line layer'
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
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Output layer')
            )
        )

    def proper_angle(self, angle: float) -> float:
        """
        Returns angles only in the range [0, 180] degrees.
        For angles bigger than 180 degrees, gives a value
        of the complementary angle to 360 degrees.
        Parameters:
        ==========
        :param angle: the input angle in radians
        Returns:
        ==========
        :return: angle in radians between 0 and pi
        """
        return angle if 0 <= angle <= pi else 2 * pi - angle

 
    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        (sink, dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT, context,
            source.fields(),
            source.wkbType(),
            source.sourceCrs()
        )

        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))

        total = 100.0 / source.featureCount() if source.featureCount() else 0
        features = source.getFeatures()

        def analyze_curve_pattern(points, start_idx, end_idx):
            """
            Analyzes curve pattern using arc length vs chord length ratios
            Returns apex point index and curve characteristics
            """
            ratios = []
            deltas = []  # Will store rate of change in ratios
            
            # Calculate arc-to-chord ratios for each point
            for i in range(start_idx, end_idx + 1):
                # Arc length (cumulative distance along points)
                arc_length = sum(
                    points[j].distance(points[j + 1])
                    for j in range(start_idx, i)
                )
                
                # Chord length (direct distance)
                chord_length = points[start_idx].distance(points[i])
                
                ratio = arc_length / chord_length if chord_length > 0 else 0
                ratios.append((i, ratio))
            
            # Calculate rate of change between consecutive ratios
            for i in range(1, len(ratios)):
                delta = ratios[i][1] - ratios[i-1][1]
                deltas.append((ratios[i][0], delta))
            
            # Look for pattern characteristic of a curve:
            # 1. Increasing positive deltas (curve entry)
            # 2. Peak (apex)
            # 3. Decreasing negative deltas (curve exit)
            
            max_delta_idx = None
            max_delta = 0
            inflection_points = []
            
            for i in range(1, len(deltas)):
                curr_delta = deltas[i][1]
                prev_delta = deltas[i-1][1]
                
                # Find where rate of change switches from increasing to decreasing
                if prev_delta > 0 and curr_delta < 0:
                    inflection_points.append(deltas[i][0])
                
                # Track maximum rate of change
                if abs(curr_delta) > max_delta:
                    max_delta = abs(curr_delta)
                    max_delta_idx = deltas[i][0]
            
            # Analyze the pattern
            if inflection_points:
                # The apex is likely near the inflection point with highest ratio
                apex_candidates = [
                    (idx, ratios[i][1])
                    for i, (idx, _) in enumerate(ratios)
                    if idx in inflection_points
                ]
                
                if apex_candidates:
                    apex_idx = max(apex_candidates, key=lambda x: x[1])[0]
                    
                    # Calculate curve characteristics
                    curve_info = {
                        'apex_idx': apex_idx,
                        'max_ratio': max(r for _, r in ratios),
                        'curve_length': arc_length,
                        'curve_deviation': max_delta,
                        'symmetry': abs(apex_idx - (start_idx + end_idx)/2) / (end_idx - start_idx)
                    }
                    
                    return curve_info
            
            return None

        for current, feature in enumerate(features):
            if feedback.isCanceled():
                break

            geom = feature.geometry()
            if not geom:
                continue

            # Get points from geometry
            points = []
            for vertex in geom.vertices():
                points.append(QgsPointXY(vertex.x(), vertex.y()))

            # Need at least 3 points to find curves
            if len(points) < 3:
                continue

            # 1. First find direction change points (as before)
            turn_points = []
            prev_angle_change = 0
            
            for i in range(1, len(points) - 1):
                azimuth1 = QgsGeometryUtils.lineAngle(
                    points[i-1].x(), points[i-1].y(),
                    points[i].x(), points[i].y()
                )
                azimuth2 = QgsGeometryUtils.lineAngle(
                    points[i].x(), points[i].y(),
                    points[i+1].x(), points[i+1].y()
                )
                
                angle_change = self.proper_angle(azimuth2 - azimuth1)
                
                if i > 1 and (
                    (prev_angle_change > 0 and angle_change < 0) or 
                    (prev_angle_change < 0 and angle_change > 0)
                ):
                    turn_points.append(i)
                
                prev_angle_change = angle_change

            # 2. Analyze sections between turn points for curve patterns
            curves = []
            for i in range(len(turn_points) - 1):
                start_idx = turn_points[i]
                end_idx = turn_points[i + 1]
                
                curve_info = analyze_curve_pattern(points, start_idx, end_idx)
                if curve_info:
                    curves.append(curve_info)

            # 3. Create output features marking curves and their characteristics
            for curve in curves:
                out_feat = QgsFeature(feature.fields())
                out_feat.setAttributes(feature.attributes())
                
                # Create point geometry for apex
                apex_point = points[curve['apex_idx']]
                out_feat.setGeometry(QgsGeometry.fromPointXY(apex_point))
                
                (sink, dest_id) = self.parameterAsSink(
                    parameters, self.OUTPUT, context,
                    source.fields(),
                    source.wkbType(),  # This copies the input geometry type (LineString)
                    source.sourceCrs()
                )

            feedback.setProgress(int(current * total))

        return {self.OUTPUT: dest_id} 