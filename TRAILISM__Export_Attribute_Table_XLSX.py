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

import os
import re

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsFeature,
    QgsFeatureRequest,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterFile,
    QgsProcessingParameterVectorLayer,
    QgsSettings,
    QgsVectorFileWriter,
    QgsVectorLayer,
)

SETTINGS_KEY_LAST_DIR = 'TRAILISM/ExportAttributeTableXlsx/last_dir'


def _safe_filename(name: str) -> str:
    name = (name or 'layer').strip()
    name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', name)
    name = name.rstrip(' .')
    return name or 'layer'


def _layer_source_dir(layer) -> str:
    """Directory of a file-based layer datasource, else ''."""
    if layer is None:
        return ''
    src = (layer.source() or '').split('|', 1)[0].strip()
    if not src or '://' in src:
        return ''
    path = os.path.normpath(src)
    if os.path.isfile(path):
        return os.path.dirname(path)
    return ''


def _default_output_dir(layer=None) -> str:
    """Prefer last-used folder, then layer source folder."""
    last = QgsSettings().value(SETTINGS_KEY_LAST_DIR, '', type=str)
    if last and os.path.isdir(last):
        return last
    if layer is not None:
        return _layer_source_dir(layer)
    from qgis.utils import iface
    if iface is None:
        return ''
    return _layer_source_dir(iface.activeLayer())


class ExportAttributeTableXlsxAlgorithm(QgsProcessingAlgorithm):
    """Export attribute table to .xlsx (data columns only, no geometry/metadata)."""

    INPUT = 'INPUT'
    OUTPUT_DIR = 'OUTPUT_DIR'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ExportAttributeTableXlsxAlgorithm()

    def name(self):
        return 'export_attribute_table_xlsx'

    def displayName(self):
        return self.tr('TRAILISM: Export Attribute Table to XLSX')

    def group(self):
        return self.tr('TRAILISM')

    def groupId(self):
        return 'trailism'

    def shortHelpString(self):
        return self.tr(
            'TRAILISM ALGORITHM:\n\n'
            'Exports the attribute table of a vector layer to Excel (.xlsx).\n'
            'Only data columns — no geometry, no WKT, no layer metadata.\n\n'
            'Minimal input:\n'
            '- Layer (defaults to the active layer when available)\n'
            '- Output folder: last-used folder first, else layer source folder\n\n'
            'File name is always the layer name + .xlsx.'
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT,
                self.tr('Input layer'),
                types=[QgsProcessing.TypeVector],
            )
        )

        default_dir = _default_output_dir() or None
        self.addParameter(
            QgsProcessingParameterFile(
                self.OUTPUT_DIR,
                self.tr('Output folder'),
                behavior=QgsProcessingParameterFile.Folder,
                optional=True,
                defaultValue=default_dir,
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        layer = self.parameterAsVectorLayer(parameters, self.INPUT, context)
        if layer is None or not layer.isValid():
            raise QgsProcessingException(self.tr('Invalid input layer'))

        out_dir = (self.parameterAsString(parameters, self.OUTPUT_DIR, context) or '').strip()
        if not out_dir:
            out_dir = _default_output_dir(layer)
        if not out_dir:
            raise QgsProcessingException(
                self.tr('Output folder is required (no last-used or layer source folder)')
            )

        out_dir = os.path.normpath(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        QgsSettings().setValue(SETTINGS_KEY_LAST_DIR, out_dir)

        out_path = os.path.join(out_dir, f'{_safe_filename(layer.name())}.xlsx')

        mem = QgsVectorLayer('None', 'attrs', 'memory')
        if not mem.isValid():
            raise QgsProcessingException(self.tr('Failed to create temporary attribute layer'))

        provider = mem.dataProvider()
        provider.addAttributes(layer.fields().toList())
        mem.updateFields()

        total = layer.featureCount()
        feats = []
        request = QgsFeatureRequest().setFlags(QgsFeatureRequest.NoGeometry)
        for i, src in enumerate(layer.getFeatures(request)):
            if feedback.isCanceled():
                return {}
            if total > 0 and i % 500 == 0:
                feedback.setProgress(int(50 * i / total))
            dst = QgsFeature(mem.fields())
            dst.setAttributes(src.attributes())
            feats.append(dst)

        ok, _ = provider.addFeatures(feats)
        if not ok:
            raise QgsProcessingException(self.tr('Failed to copy attributes'))

        if feedback.isCanceled():
            return {}

        options = QgsVectorFileWriter.SaveVectorOptions()
        options.driverName = 'XLSX'
        options.fileEncoding = 'UTF-8'
        options.layerOptions = ['HEADERS=FORCE']

        err_code, err_msg, written_path, _layer_name = QgsVectorFileWriter.writeAsVectorFormatV3(
            mem,
            out_path,
            context.transformContext(),
            options,
        )
        if err_code != QgsVectorFileWriter.NoError:
            raise QgsProcessingException(
                self.tr(f'XLSX write failed ({err_code}): {err_msg}')
            )

        feedback.setProgress(100)
        feedback.pushInfo(self.tr(f'Saved: {written_path}'))
        return {self.OUTPUT: written_path}
