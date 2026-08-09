"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
Raster analogue of tree-top (local maxima) detection from:
Tiede, D., Hochleitner, G., Blaschke, T. (2005). A fully GIS-based workflow for tree
identification and tree crown delineation using laser scanning.
CD = a + b*TH (defaults a=1.54, b=0.123); search radius R = CD/2 at each focal cell.
Peak localization only — no crown delineation.
"""

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterVectorDestination,
    QgsProcessingParameterExtent,
    QgsProcessingException,
    QgsPointXY,
    QgsVectorFileWriter,
    QgsFields,
    QgsField,
    QgsFeature,
    QgsGeometry,
    QgsWkbTypes,
)
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
import math
import os
import tempfile

import numpy as np

try:
    from osgeo import gdal
    GDAL_AVAILABLE = True
except ImportError:
    GDAL_AVAILABLE = False


def _disk_footprint(radius_px):
    """Boolean 2D footprint: disk of integer radius radius_px."""
    r = int(max(0, radius_px))
    yy, xx = np.ogrid[-r : r + 1, -r : r + 1]
    return (xx * xx + yy * yy) <= r * r


def _neighbors_in_disk_meters(
    height, neighbor_valid, r0, c0, r_m, cell_size_x, cell_size_y
):
    """Local maximum within disk: no neighbor strictly higher; ties broken by row-major order (r,c)."""
    h0 = height[r0, c0]
    if not np.isfinite(h0):
        return False
    max_px = int(math.ceil(r_m / min(cell_size_x, cell_size_y))) + 1
    h, w = height.shape
    for dr in range(-max_px, max_px + 1):
        for dc in range(-max_px, max_px + 1):
            if dr == 0 and dc == 0:
                continue
            rr, cc = r0 + dr, c0 + dc
            if rr < 0 or rr >= h or cc < 0 or cc >= w:
                continue
            dist_m = math.hypot(dr * cell_size_x, dc * cell_size_y)
            if dist_m > r_m + 1e-9:
                continue
            if not neighbor_valid[rr, cc]:
                continue
            hv = height[rr, cc]
            if hv > h0:
                return False
            if hv == h0 and (rr, cc) < (r0, c0):
                return False
    return True


def _find_seeds_tiede(
    height,
    core_valid,
    neighbor_valid,
    a,
    b,
    max_crown_diameter_m,
    cell_size_x,
    cell_size_y,
):
    """
    Per-cell: CD = min(a + b*TH, max_crown_diameter_m), R_m = CD/2 (core_valid cells only).
    Coarse candidates: focal max on disk of max R (pixels).
    Peaks: local maximum in disk of radius R_m (meters); ties broken by row-major (r,c) first.
    Returns: seeds_rc, r_m_map.
    """
    th = np.where(core_valid, height, np.nan)
    cd = a + b * np.maximum(th, 0.0)
    cd = np.minimum(cd, float(max_crown_diameter_m))
    r_m_map = np.where(core_valid, cd * 0.5, 0.0)

    max_r_m = float(max_crown_diameter_m) * 0.5
    max_px = int(math.ceil(max_r_m / min(cell_size_x, cell_size_y)))
    if max_px < 1:
        max_px = 1

    try:
        from scipy import ndimage
    except Exception as e:
        raise QgsProcessingException(
            "scipy.ndimage is required for seed detection. Install scipy in the QGIS Python, "
            "or if import fails with NumPy ufunc errors: python -m pip install -U numpy scipy"
        ) from e

    footprint = _disk_footprint(max_px)
    h_for_max = np.where(core_valid, height, -np.inf)
    local_max_vals = ndimage.maximum_filter(
        h_for_max, footprint=footprint, mode="constant", cval=-np.inf
    )
    candidate_mask = core_valid & np.isfinite(height) & (
        height >= local_max_vals - 1e-9
    )

    seeds_rc = []
    h, w = height.shape
    for r0 in range(h):
        for c0 in range(w):
            if not candidate_mask[r0, c0]:
                continue
            rm = r_m_map[r0, c0]
            if rm <= 0:
                continue
            if _neighbors_in_disk_meters(
                height, neighbor_valid, r0, c0, rm, cell_size_x, cell_size_y
            ):
                seeds_rc.append((r0, c0))

    return seeds_rc, r_m_map


class TreePeaksTiedeStyleAlgorithm(QgsProcessingAlgorithm):
    INPUT_RASTER = "INPUT_RASTER"
    INPUT_DTM = "INPUT_DTM"
    EXTENT = "EXTENT"
    COEFF_A = "COEFF_A"
    COEFF_B = "COEFF_B"
    MAX_CROWN_DIAMETER_M = "MAX_CROWN_DIAMETER_M"
    MIN_TREE_HEIGHT_M = "MIN_TREE_HEIGHT_M"
    EDGE_PAD_M = "EDGE_PAD_M"
    OUTPUT_PEAKS = "OUTPUT_PEAKS"

    def tr(self, string):
        return QCoreApplication.translate("Processing", string)

    def createInstance(self):
        return TreePeaksTiedeStyleAlgorithm()

    def name(self):
        return "tree_peaks_tiede_style"

    def displayName(self):
        return "TRAILISM: Tree peaks (Tiede local maxima)"

    def shortHelpString(self):
        return (
            "Tree top localization using Tiede-style local maxima: CD = min(a + b*TH, cap), "
            "search radius R = CD/2 at each focal cell (coarse max filter, then per-cell disk; ties: row-major).\n\n"
            "Edge padding (meters) extends the grid with low synthetic cells so disks near the "
            "raster border are fully sampled; outputs are only written for the original extent.\n\n"
            "No crown delineation — peak points only. Defaults a=1.54, b=0.123 from Tiede et al. (2005).\n"
            "Input: DSM or CHM/nDOM; optional DTM for nDOM = DSM − DTM; optional extent.\n"
            "Requires scipy.ndimage for the coarse max filter."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_RASTER,
                self.tr("Input raster (DSM or canopy height CHM/nDOM)"),
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DTM,
                self.tr("DTM (optional) — if set, nDOM = DSM − DTM"),
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterExtent(
                self.EXTENT,
                self.tr("Extent (optional) — limit processing to this area"),
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.COEFF_A,
                self.tr("Crown model: a in CD = a + b·TH (meters)"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.54,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.COEFF_B,
                self.tr("Crown model: b in CD = a + b·TH"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.123,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_CROWN_DIAMETER_M,
                self.tr("Maximum crown diameter cap (meters)"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=10.0,
                minValue=0.1,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MIN_TREE_HEIGHT_M,
                self.tr("Minimum canopy height — preliminary mask (meters)"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=5.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.EDGE_PAD_M,
                self.tr("Edge padding (meters) — pad grid for full disk sampling"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=10.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterVectorDestination(
                self.OUTPUT_PEAKS,
                self.tr("Output tree peaks (local maxima, row-major tie-break)"),
                type=QgsProcessing.TypeVectorPoint,
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        if not GDAL_AVAILABLE:
            raise QgsProcessingException("GDAL Python bindings are required.")

        raster_layer = self.parameterAsRasterLayer(parameters, self.INPUT_RASTER, context)
        if raster_layer is None or not raster_layer.isValid():
            raise QgsProcessingException("Invalid input raster.")

        dtm_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DTM, context)
        a = self.parameterAsDouble(parameters, self.COEFF_A, context)
        b = self.parameterAsDouble(parameters, self.COEFF_B, context)
        max_cd = self.parameterAsDouble(parameters, self.MAX_CROWN_DIAMETER_M, context)
        min_h = self.parameterAsDouble(parameters, self.MIN_TREE_HEIGHT_M, context)
        pad_m = self.parameterAsDouble(parameters, self.EDGE_PAD_M, context)

        peaks_path = self.parameterAsOutputLayer(parameters, self.OUTPUT_PEAKS, context)
        if not peaks_path:
            raise QgsProcessingException("Output tree peaks path is required.")

        crs = raster_layer.crs()
        raster_extent = raster_layer.extent()
        extent_rect = self.parameterAsExtent(parameters, self.EXTENT, context, crs)
        if extent_rect and not extent_rect.isEmpty() and not raster_extent.intersects(extent_rect):
            raise QgsProcessingException("Provided extent does not overlap the input raster extent.")

        temp_dir = tempfile.mkdtemp()
        temp_files = []

        try:
            input_source = raster_layer.source()
            if dtm_layer is not None:
                feedback.setProgressText("Computing nDOM (DSM − DTM)...")
                ndom_path = os.path.join(temp_dir, "nDOM.tif")
                temp_files.append(ndom_path)
                e_dsm = QgsRasterCalculatorEntry()
                e_dsm.ref = "dsm@1"
                e_dsm.raster = raster_layer
                e_dsm.bandNumber = 1
                e_dtm = QgsRasterCalculatorEntry()
                e_dtm.ref = "dtm@1"
                e_dtm.raster = dtm_layer
                e_dtm.bandNumber = 1
                calc = QgsRasterCalculator(
                    "dsm@1 - dtm@1",
                    ndom_path,
                    "GTiff",
                    raster_layer.extent(),
                    raster_layer.crs(),
                    raster_layer.width(),
                    raster_layer.height(),
                    [e_dsm, e_dtm],
                    context.transformContext(),
                )
                res = calc.processCalculation(feedback)
                if res != QgsRasterCalculator.Success:
                    raise QgsProcessingException(
                        f"nDOM failed: {calc.lastError() or 'unknown'}"
                    )
                input_source = ndom_path

            if extent_rect and not extent_rect.isEmpty():
                feedback.setProgressText("Clipping to extent...")
                clipped_path = os.path.join(temp_dir, "clipped.tif")
                temp_files.append(clipped_path)
                x_min = extent_rect.xMinimum()
                x_max = extent_rect.xMaximum()
                y_min = extent_rect.yMinimum()
                y_max = extent_rect.yMaximum()
                proj_win = [x_min, y_max, x_max, y_min]
                translate_kwargs = {"format": "GTiff", "projWin": proj_win}
                auth_id = crs.authid() if crs.isValid() else ""
                if auth_id:
                    translate_kwargs["projWinSRS"] = auth_id
                ds_clip = gdal.Translate(
                    clipped_path,
                    input_source,
                    options=gdal.TranslateOptions(**translate_kwargs),
                )
                if ds_clip is None:
                    raise QgsProcessingException(
                        "Extent clipping failed (GDAL Translate / projWin)."
                    )
                ds_clip = None
                input_source = clipped_path

            ds = gdal.Open(input_source)
            if ds is None:
                raise QgsProcessingException("Could not open height raster.")
            band = ds.GetRasterBand(1)
            height_raw = band.ReadAsArray().astype(np.float64)
            nodata = band.GetNoDataValue()
            gt = ds.GetGeoTransform()
            ds = None

            xmin = gt[0]
            ymax = gt[3]
            cell_size_x = abs(gt[1])
            cell_size_y = abs(gt[5])

            if nodata is None:
                valid_finite = np.isfinite(height_raw)
            else:
                valid_finite = np.isfinite(height_raw) & (height_raw != nodata)

            low_fill = min_h - 1.0
            height_orig = np.where(valid_finite, height_raw, low_fill)
            core_candidate = valid_finite & (height_raw >= min_h)

            h0, w0 = height_orig.shape
            pad_px = max(0, int(math.ceil(pad_m / cell_size_x)))
            pad_py = max(0, int(math.ceil(pad_m / cell_size_y)))

            height = np.pad(
                height_orig,
                ((pad_py, pad_py), (pad_px, pad_px)),
                mode="constant",
                constant_values=low_fill,
            )
            core_valid = np.pad(
                core_candidate,
                ((pad_py, pad_py), (pad_px, pad_px)),
                mode="constant",
                constant_values=False,
            )
            neighbor_valid = np.ones_like(height, dtype=bool)

            feedback.setProgressText("Detecting coarse candidates and peaks (Tiede R, row-major ties)...")
            seeds_rc, r_m_map = _find_seeds_tiede(
                height,
                core_valid,
                neighbor_valid,
                a,
                b,
                max_cd,
                cell_size_x,
                cell_size_y,
            )

            def _in_core(r, c):
                return (
                    pad_py <= r < pad_py + h0
                    and pad_px <= c < pad_px + w0
                )

            seeds_rc = [(r, c) for r, c in seeds_rc if _in_core(r, c)]

            feedback.setProgressText("Writing tree peak points...")
            fld_peak = QgsFields()
            fld_peak.append(QgsField("peak_id", QVariant.Int))
            fld_peak.append(QgsField("tree_height", QVariant.Double))
            fld_peak.append(QgsField("cd_eff_m", QVariant.Double))
            fld_peak.append(QgsField("r_search_m", QVariant.Double))
            drv_peak = (
                "ESRI Shapefile"
                if peaks_path.lower().endswith(".shp")
                else "GPKG"
            )
            writer_p = QgsVectorFileWriter(
                peaks_path, "UTF-8", fld_peak, QgsWkbTypes.Point, crs, drv_peak
            )
            if writer_p.hasError():
                raise QgsProcessingException(writer_p.errorMessage())
            for k, (r0, c0) in enumerate(seeds_rc):
                c_orig = c0 - pad_px
                r_orig = r0 - pad_py
                xs = xmin + (c_orig + 0.5) * cell_size_x
                ys = ymax - (r_orig + 0.5) * cell_size_y
                th = float(height_raw[r_orig, c_orig])
                cd_eff = float(min(a + b * th, max_cd))
                r_s = float(r_m_map[r0, c0])
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(xs, ys)))
                feat.setAttributes([k + 1, th, cd_eff, r_s])
                writer_p.addFeature(feat)
            del writer_p

            result = {self.OUTPUT_PEAKS: peaks_path}
            feedback.pushInfo(f"Wrote {len(seeds_rc)} tree peaks to {peaks_path}")

            feedback.setProgress(100)
            return result
        finally:
            for f in temp_files:
                try:
                    if os.path.isfile(f):
                        os.remove(f)
                except OSError:
                    pass
            try:
                os.rmdir(temp_dir)
            except OSError:
                pass
