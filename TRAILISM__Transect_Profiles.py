# pyright: reportMissingImports=false
from qgis.PyQt.QtCore import QVariant, QUrl  # type: ignore[import]
from qgis.PyQt.QtGui import QDesktopServices  # type: ignore[import]
import math
import os
from qgis.core import (  # type: ignore[import]
	QgsProcessing,
	QgsProcessingAlgorithm,
	QgsProcessingParameterFeatureSource,
	QgsProcessingParameterNumber,
	QgsProcessingParameterBoolean,
	QgsProcessingParameterEnum,
	QgsProcessingParameterString,
	QgsProcessingParameterField,
	QgsProcessingParameterFeatureSink,
	QgsProcessingParameterRasterLayer,
	QgsProcessingException,
	QgsFields,
	QgsField,
	QgsFeature,
	QgsFeatureSink,
	QgsGeometry,
	QgsPointXY,
	QgsWkbTypes,
	QgsCoordinateTransform,
	QgsProject,
	QgsRaster,
	QgsMarkerSymbol,
	QgsLineSymbol,
	QgsSingleSymbolRenderer,
	QgsRuleBasedRenderer,
)


# Unified matplotlib bbox style for in-plot labels
MPL_BBOX = dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.5)

# Safe filename token from arbitrary text
def _sanitize_filename_component(txt):
	try:
		s = str(txt)
	except Exception:
		return "feat"
	safe = []
	for ch in s:
		if ch.isalnum() or ch in ('_', '-', '.'):
			safe.append(ch)
		elif ch.isspace():
			safe.append('_')
		else:
			safe.append('_')
	out = ''.join(safe).strip('_')
	return out if out else "feat"

# This script defines a QGIS Processing tool that samples cross-sections
# perpendicular to a trail/centerline and compares a simple "design" profile
# to a DEM (Digital Elevation Model). From that, it finds where side slopes
# meet natural ground, computes cut/fill areas, and optionally draws an SVG
# profile diagram for quick QA.
#
# Key terms in plain words:
# - DEM: a raster grid where each pixel stores ground elevation.
# - Station: distance along the input line (0 at the start of each part).
# - Normal: the direction perpendicular to the line's local tangent.
# - +s coordinate: distance along the normal. Here, +s is defined to point to the
#   geometric RIGHT side of the line; -s points to the LEFT side.
# - Trail edges: s = -W/2 (left) and s = +W/2 (right), with W = trail width.
# - beta_1: side slope angle (degrees) used to extend from the trail edge to meet
#   existing ground. tan(beta_1) is the slope ratio (rise over run).
# - C point: the intersection of the design side slope with the DEM profile
#   on that side (toe/heel of slope). We compute one C per side if it exists.
# - cut/fill: whether the design sits below (cut) or above (fill) the natural ground.

class GenerateProfilesAndSlopeIntersections(QgsProcessingAlgorithm):
	LINE_LAYER = "LINE_LAYER"
	DEM = "DEM"
	STEP_M = "STEP_M"
	WIDTH_FIELD = "WIDTH_FIELD"
	WIDTH_FALLBACK = "WIDTH_FALLBACK"
	BETA1_DEG = "BETA1_DEG"
	PROFILE_RANGE_M = "PROFILE_RANGE_M"
	AXIS_DZ = "AXIS_DZ"
	AXIS_FULL_CUT = "AXIS_FULL_CUT"
	AXIS_FULL_FILL = "AXIS_FULL_FILL"
	DEM_OFFSET = "DEM_OFFSET"
	SURF_HEIGHT = "SURF_HEIGHT"
	SURF_CROSSFALL_DEG = "SURF_CROSSFALL_DEG"
	SURF_TOP_TILT_DEG = "SURF_TOP_TILT_DEG"
	SURF_TOP_CROSSFALL_MODE = "SURF_TOP_CROSSFALL_MODE"
	SURF_HEIGHT_FIELD = "SURF_HEIGHT_FIELD"
	SVG_PX_PER_M = "SVG_PX_PER_M"

	OUT_POINTS = "OUT_POINTS"
	OUT_BC_SEGMENTS = "OUT_BC_SEGMENTS"
	OUT_AREAS = "OUT_AREAS"
	OUT_CHAIN_C = "OUT_CHAIN_C"
	OUT_CHAIN_EDGES = "OUT_CHAIN_EDGES"  # kept for compatibility; returns same as OUT_CHAIN_C
	OUT_TRANSECTS = "OUT_TRANSECTS"      # kept for compatibility; returns same as OUT_BC_SEGMENTS

	# SVG preview controls
	SVG_OPEN = "SVG_OPEN"
	SVG_STATIONS = "SVG_STATIONS"
	SVG_ALL_STATIONS = "SVG_ALL_STATIONS"

	def initAlgorithm(self, config=None):
		# Define all input parameters shown in the QGIS Processing panel.
		# QGIS will validate and pass these to processAlgorithm.
		self.addParameter(QgsProcessingParameterFeatureSource(
			self.LINE_LAYER, "Input trail centerline (line)", [QgsProcessing.TypeVectorLine]
		))
		self.addParameter(QgsProcessingParameterRasterLayer(self.DEM, "DEM raster"))
		self.addParameter(
			QgsProcessingParameterNumber(
				self.BETA1_DEG,
				"Side slope angle beta_1 (degrees; 0–89.999). tan(beta_1) gives rise/run; legend shows nearest simple ratio.",
				QgsProcessingParameterNumber.Double,
				defaultValue=33.7, # 2:3 slope ratio
				minValue=0.0,
				maxValue=89.999,
			)
		)
		self.addParameter(
			QgsProcessingParameterNumber(
				self.STEP_M,
				"Station spacing along axis (meters; ≥ 1). Stations are multiples of this step.",
				QgsProcessingParameterNumber.Double,
				defaultValue=10.0,
				minValue=1.0,
			)
		)
		self.addParameter(
			QgsProcessingParameterField(
				self.WIDTH_FIELD,
				"Width attribute field (meters). If set and valid per feature, overrides the default width.",
				parentLayerParameterName=self.LINE_LAYER,
				type=QgsProcessingParameterField.Numeric,
				optional=True,
			)
		)
		self.addParameter(
			QgsProcessingParameterNumber(
				self.WIDTH_FALLBACK,
				"Default trail width (meters; ≥ 0). Used when width field is missing/invalid.",
				QgsProcessingParameterNumber.Double,
				defaultValue=2.0,
				minValue=0.0,
			)
		)
		self.addParameter(
			QgsProcessingParameterNumber(
				self.PROFILE_RANGE_M,
				"Profile half-range R (meters; ≥ 1). DEM is sampled from -R (left) to +R (right) around axis.",
				QgsProcessingParameterNumber.Double,
				defaultValue=2.5,
				minValue=1.0,
			)
		)
		# Vertical offset applied to DEM everywhere (computations and plots)
		self.addParameter(
			QgsProcessingParameterNumber(
				self.DEM_OFFSET,
				"DEM vertical offset (meters; -1 to 0). Applied globally to DEM for all calculations and plots.",
				QgsProcessingParameterNumber.Double,
				defaultValue=-0.15,
				minValue=-1,
				maxValue=0,
			)
		)
		# Optional axis elevation offset: z0 := DEM + dz
		self.addParameter(
			QgsProcessingParameterNumber(
				self.AXIS_DZ,
				"Axis elevation offset dz (meters; -2.5 to 2.5). z0 := DEM + dz unless overridden by full cut/fill mode.",
				QgsProcessingParameterNumber.Double,
				defaultValue=-0.1,
				minValue=-2.5,
				maxValue=2.5,
			)
		)
		# Axis z override modes (ignore dz when enabled)
		self.addParameter(QgsProcessingParameterBoolean(self.AXIS_FULL_CUT, "Axis mode: full cut — set z0 to the lower DEM edge (min). Ignores dz.", defaultValue=False))
		self.addParameter(QgsProcessingParameterBoolean(self.AXIS_FULL_FILL, "Axis mode: full fill — set z0 to the higher DEM edge (max). Ignores dz.", defaultValue=False))

		# SVG surface overlay parameters (render-only): thickness and crossfall
		self.addParameter(
			QgsProcessingParameterNumber(
				self.SURF_HEIGHT,
				"Surface thickness height h_A (meters; 0–1). Visual overlay only; does not affect cut/fill.",
				QgsProcessingParameterNumber.Double,
				defaultValue=0.2,
				minValue=0.0,
				maxValue=1.0,
			)
		)

		self.addParameter(
			QgsProcessingParameterField(
				self.SURF_HEIGHT_FIELD,
				"Surface thickness field (meters). If present per feature, overrides h_A.",
				parentLayerParameterName=self.LINE_LAYER,
				type=QgsProcessingParameterField.Numeric,
				optional=True,
			)
		)

		self.addParameter(
			QgsProcessingParameterNumber(
				self.SURF_CROSSFALL_DEG,
				"Side face angle to horizontal (degrees; 20–90). 90° = vertical sides (rectangle). Visual overlay only.",
				QgsProcessingParameterNumber.Double,
				defaultValue=45.0,
				minValue=20.0,
				maxValue=90.0,
			)
		)

		# Visual-only: surface top crossfall mode (does not affect calculations)
		self.addParameter(
			QgsProcessingParameterEnum(
				self.SURF_TOP_CROSSFALL_MODE,
				"Surface top crossfall (visual only)",
				options=["Flat", "Inward (to uphill)", "Outward (to downhill)"],
				allowMultiple=False,
				defaultValue=2,
			)
		)
		# Visual-only: crossfall tilt angle degrees when Inward/Outward is chosen
		self.addParameter(
			QgsProcessingParameterNumber(
				self.SURF_TOP_TILT_DEG,
				"Surface top crossfall tilt q (degrees; 0–7). Used only when crossfall ≠ Flat.",
				QgsProcessingParameterNumber.Integer,
				defaultValue=3,
				minValue=0,
				maxValue=7,
			)
		)

		# SVG preview (explicit stations or all stations)
		self.addParameter(QgsProcessingParameterBoolean(self.SVG_OPEN, "Open SVG after creation", defaultValue=True))
		self.addParameter(QgsProcessingParameterBoolean(self.SVG_ALL_STATIONS, "SVG: use all stations (max 100 SVGs per run)", defaultValue=False))
		self.addParameter(QgsProcessingParameterString(self.SVG_STATIONS, "SVG stations to plot (meters, comma-separated; exact match). Overall SVG cap: 100 per run.", defaultValue="", multiLine=False, optional=True))
		self.addParameter(
			QgsProcessingParameterNumber(
				self.SVG_PX_PER_M,
				"SVG scale: pixels per meter (30–300). Controls both axes; aspect is 1:1.",
				QgsProcessingParameterNumber.Double,
				defaultValue=180.0,
				minValue=10.0,
				maxValue=300.0,
			)
		)

		# Outputs
		pt_fields = QgsFields()
		pt_fields.append(QgsField("station_m", QVariant.Double))
		pt_fields.append(QgsField("ptype", QVariant.String))  # edge_left, edge_right, slope_left, slope_right
		pt_fields.append(QgsField("side", QVariant.String))   # left | right | center
		pt_fields.append(QgsField("cut_fill", QVariant.String)) # cut | fill | edge
		pt_fields.append(QgsField("slope_run", QVariant.Double)) # horizontal distance from trail edge to C
		pt_fields.append(QgsField("z", QVariant.Double))      # DEM z at point (edges use z0)
		pt_fields.append(QgsField("z_base", QVariant.Double)) # z0 at axis
		pt_fields.append(QgsField("z_diff", QVariant.Double)) # z - z0
		pt_fields.append(QgsField("svg_sel", QVariant.Int))   # 1 if selected for SVG
		# Aggregated per-station areas (set on axis_center only)
		pt_fields.append(QgsField("area_cut_total", QVariant.Double))
		pt_fields.append(QgsField("area_fill_total", QVariant.Double))
		self.addParameter(QgsProcessingParameterFeatureSink(self.OUT_POINTS, "Points: axis center, trail edges, slope toes/heels (with areas)", defaultValue=QgsProcessing.TEMPORARY_OUTPUT))

		seg_fields = QgsFields()
		seg_fields.append(QgsField("station_m", QVariant.Double))
		seg_fields.append(QgsField("side", QVariant.String))
		seg_fields.append(QgsField("cut_fill", QVariant.String))
		seg_fields.append(QgsField("slope_run", QVariant.Double))
		seg_fields.append(QgsField("dz_BC", QVariant.Double))
		# Combine BC and transect lines in one layer
		seg_fields.append(QgsField("linetype", QVariant.String))  # bc_line | transect_line
		seg_fields.append(QgsField("length_m", QVariant.Double))  # for transect_line; None for bc_line
		self.addParameter(QgsProcessingParameterFeatureSink(self.OUT_BC_SEGMENTS, "BC segments per side (line)", defaultValue=QgsProcessing.TEMPORARY_OUTPUT))

		chain_fields = QgsFields()

		# Areas per station are merged into the points output (axis_center feature)

		# Chains (lines): C_upper, C_lower, trail edges (left/right)
		chain_fields = QgsFields()
		self.addParameter(QgsProcessingParameterFeatureSink(self.OUT_CHAIN_C, "Chains: trail edges and C (line)", defaultValue=QgsProcessing.TEMPORARY_OUTPUT))

	def processAlgorithm(self, parameters, context, feedback):
		# This is the main entry point. QGIS calls this with the user's inputs.
		# We iterate each feature/part of the input line, walk stations along it,
		# build a cross-section profile from the DEM, compare to the design,
		# find intersections and areas, emit outputs, and optionally plot an SVG.
		line_source = self.parameterAsSource(parameters, self.LINE_LAYER, context)
		if line_source is None:
			raise QgsProcessingException("Invalid line layer")

		step_m = float(self.parameterAsDouble(parameters, self.STEP_M, context))
		width_field = self.parameterAsString(parameters, self.WIDTH_FIELD, context)
		width_fallback = float(self.parameterAsDouble(parameters, self.WIDTH_FALLBACK, context))
		dem_layer = self.parameterAsRasterLayer(parameters, self.DEM, context)
		if dem_layer is None:
			raise QgsProcessingException("Invalid DEM raster")
		beta1_deg = float(self.parameterAsDouble(parameters, self.BETA1_DEG, context))
		beta1_rad = math.radians(max(0.0, min(89.999, beta1_deg)))
		R = float(self.parameterAsDouble(parameters, self.PROFILE_RANGE_M, context))
		dz_axis = float(self.parameterAsDouble(parameters, self.AXIS_DZ, context))
		dem_offset = float(self.parameterAsDouble(parameters, self.DEM_OFFSET, context))
		surf_height = float(self.parameterAsDouble(parameters, self.SURF_HEIGHT, context))
		surf_top_mode = int(self.parameterAsEnum(parameters, self.SURF_TOP_CROSSFALL_MODE, context)) if hasattr(self, 'SURF_TOP_CROSSFALL_MODE') else 0
		surf_top_tilt_deg = int(self.parameterAsInt(parameters, self.SURF_TOP_TILT_DEG, context)) if hasattr(self, 'SURF_TOP_TILT_DEG') else 4
		surf_height_field = self.parameterAsString(parameters, self.SURF_HEIGHT_FIELD, context) or ""
		surf_xfall_deg = float(self.parameterAsDouble(parameters, self.SURF_CROSSFALL_DEG, context))
		# Clamp to avoid degenerate geometry and tan(90°) handling
		_sdeg = max(20.0, min(90.0, surf_xfall_deg))
		vertical_sides = (abs(_sdeg - 90.0) <= 1e-9)
		surf_xfall_rad = math.radians(_sdeg if not vertical_sides else 89.999)
		axis_full_cut = bool(self.parameterAsBool(parameters, self.AXIS_FULL_CUT, context))
		axis_full_fill = bool(self.parameterAsBool(parameters, self.AXIS_FULL_FILL, context))
		svg_open = bool(self.parameterAsBool(parameters, self.SVG_OPEN, context))
		svg_all = bool(self.parameterAsBool(parameters, self.SVG_ALL_STATIONS, context))
		svg_stations_raw = self.parameterAsString(parameters, self.SVG_STATIONS, context) or ""
		svg_px_per_m = float(self.parameterAsDouble(parameters, self.SVG_PX_PER_M, context))
		requested_stations = []
		if svg_stations_raw.strip():
			for tok in svg_stations_raw.split(','):
				try:
					val = float(tok.strip())
					requested_stations.append(val)
				except Exception:
					pass

		# DEM access
		# We prepare helpers to sample elevation values at arbitrary map points.
		MAX_SVGS = 100
		svg_count = 0
		dp = dem_layer.dataProvider()
		nodata = dp.sourceNoDataValue(1) if hasattr(dp, "sourceNoDataValue") else None
		raster_units_per_pixel_x = getattr(dem_layer.rasterUnitsPerPixelX, "__call__", None)
		px = dem_layer.rasterUnitsPerPixelX() if raster_units_per_pixel_x else None
		if step_m <= 0.0:
			h = max(0.1, float(px) if px else 1.0)
		else:
			h = max(0.05, step_m)

		# Derived steps for DEM base sampling and evaluation spacing
		dem_step = max(0.1, float(px) if px else 1.0)
		eval_step = dem_step

		# CRS transform layer -> DEM
		# If the line and the DEM are in different coordinate reference systems (CRS),
		# transform points to the DEM's CRS before sampling elevation.
		src_crs = line_source.sourceCrs()
		dem_crs = dem_layer.crs()
		xform = QgsCoordinateTransform(src_crs, dem_crs, QgsProject.instance())

		def sample_dem(pt_src_xy):
			# Sample the DEM at a single point.
			# First try a fast sample via dataProvider.sample; if it fails or returns
			# nodata, fall back to a slower identify call.
			pt = xform.transform(pt_src_xy)
			val = dp.sample(pt, 1)
			if val is not None:
				try:
					zv = float(val)
				except Exception:
					zv = None
				if zv is not None and (nodata is None or zv != nodata):
					return zv
			res = dp.identify(pt, QgsRaster.IdentifyFormatValue)
			if not res or not res.isValid():
				return None
			vals = res.results()
			z = vals.get(1)
			if z is None:
				return None
			try:
				zv2 = float(z)
			except Exception:
				return None
			if nodata is not None and zv2 == nodata:
				return None
			return zv2

		def interp_point(geom_obj, dist_along):
			# Get a coordinate at a given distance along the line. QGIS provides
			# geometry.interpolate(d) which returns a small geometry we convert to a point.
			g = geom_obj.interpolate(dist_along)
			if g is None:
				return None
			try:
				if hasattr(g, "isNull") and g.isNull():
					return None
				if g.isEmpty():
					return None
			except Exception:
				pass
			try:
				pt = g.asPoint()
				return QgsPointXY(pt.x(), pt.y())
			except Exception:
				return None

		# (Removed name sanitization; filenames will use feature id)

		# Prepare sinks
		# Sinks are the output layers we will fill: points, line segments, areas, and chains.
		pt_fields = QgsFields()
		pt_fields.append(QgsField("station_m", QVariant.Double))
		pt_fields.append(QgsField("ptype", QVariant.String))
		pt_fields.append(QgsField("side", QVariant.String))
		pt_fields.append(QgsField("cut_fill", QVariant.String))
		pt_fields.append(QgsField("slope_run", QVariant.Double))
		pt_fields.append(QgsField("z", QVariant.Double))
		pt_fields.append(QgsField("z_base", QVariant.Double))
		pt_fields.append(QgsField("z_diff", QVariant.Double))
		pt_fields.append(QgsField("svg_sel", QVariant.Int))
		# Aggregated per-station areas (set on axis_center only)
		pt_fields.append(QgsField("area_cut_total", QVariant.Double))
		pt_fields.append(QgsField("area_fill_total", QVariant.Double))
		points_sink, points_id = self.parameterAsSink(
			parameters, self.OUT_POINTS, context, pt_fields, QgsWkbTypes.Point, line_source.sourceCrs()
		)
		points_added = 0

		# (Removed experimental canvas auto-load helper)

		seg_fields = QgsFields()
		seg_fields.append(QgsField("station_m", QVariant.Double))
		seg_fields.append(QgsField("side", QVariant.String))
		seg_fields.append(QgsField("cut_fill", QVariant.String))
		seg_fields.append(QgsField("slope_run", QVariant.Double))
		seg_fields.append(QgsField("dz_BC", QVariant.Double))
		# Ensure merged schema matches BC + transect lines
		seg_fields.append(QgsField("linetype", QVariant.String))  # 'bc_line' | 'transect_line'
		seg_fields.append(QgsField("length_m", QVariant.Double))  # only for transect_line
		seg_fields.append(QgsField("name", QVariant.String))      # optional source feature name
		segments_sink, segments_id = self.parameterAsSink(
			parameters, self.OUT_BC_SEGMENTS, context, seg_fields, QgsWkbTypes.LineString, line_source.sourceCrs()
		)

		# No profile sink (profiles are not emitted)
		# We only keep the per-sample arrays in memory for computation and plotting.

		# No separate areas sink anymore; areas go into points features (axis_center)

		# Chain sinks
		# Chains are polylines that connect C points (upper/lower) and the trail edges across stations.
		chain_fields2 = QgsFields()
		# Store which side the chain represents: 'left' or 'right'
		chain_fields2.append(QgsField("side", QVariant.String))
		# ctype distinguishes chain kind: 'trail_edge' or 'slope_edge'
		chain_fields2.append(QgsField("ctype", QVariant.String))
		chain_c_sink, chain_c_id = self.parameterAsSink(
			parameters, self.OUT_CHAIN_C, context, chain_fields2, QgsWkbTypes.LineString, line_source.sourceCrs()
		)
		# Use the same sink for edges and C chains so both go into the same output layer
		chain_edges_sink, chain_edges_id = chain_c_sink, chain_c_id
		# Remove separate transect sink; we will emit transect lines into segments_sink
		transect_sink, transect_id = segments_sink, segments_id

		feature_count = line_source.featureCount() or 0
		for idx, feat in enumerate(line_source.getFeatures()):
			if feedback.isCanceled():
				break
			if feature_count:
				feedback.setProgress(int(100.0 * idx / feature_count))

			geom = feat.geometry()
			if not geom or geom.isEmpty():
				continue

			parts = []
			if geom.isMultipart():
				for poly in geom.asMultiPolyline():
					parts.append(poly)
			else:
				parts.append(geom.asPolyline())

			# Determine longest part index for this feature
			part_lengths = []
			for p in parts:
				try:
					pl = QgsGeometry.fromPolylineXY([QgsPointXY(pt.x(), pt.y()) for pt in p]).length()
					part_lengths.append(pl)
				except Exception:
					part_lengths.append(0.0)
			longest_part_index = max(range(len(parts)), key=lambda i: part_lengths[i]) if parts else 0
			svg_done_for_feature = False
			matched_stations_for_feature = set()
			# Plot/file base name: try 'name'/'Name' field; fallback to feature id
			_plot_name = None
			try:
				_field_names = list(feat.fields().names())
				_name_field = next((n for n in _field_names if isinstance(n, str) and n.lower() == 'name'), None)
				if _name_field is not None:
					_val = feat[_name_field]
					if _val is not None and str(_val).strip():
						_plot_name = _sanitize_filename_component(_val)
			except Exception:
				pass
			if not _plot_name:
				_plot_name = f"feat_{feat.id()}"

			# width for Feature (W)
			# If a field is provided, use it; otherwise fall back to the default width.
			W = width_fallback
			if width_field and width_field in feat.fields().names():
				val = feat[width_field]
				if val is not None:
					try:
						W = float(val)
					except Exception:
						W = width_fallback
			W = max(0.0, W)

			# Per-feature surface thickness override from field, else parameter
			surf_height_eff = surf_height
			if surf_height_field and surf_height_field in feat.fields().names():
				try:
					_valh = feat[surf_height_field]
					if _valh is not None:
						sh = float(_valh)
						if sh >= 0.0:
							surf_height_eff = sh
				except Exception:
					pass

			for part_idx, part in enumerate(parts):
				if not part or len(part) < 2:
					continue
				part_geom = QgsGeometry.fromPolylineXY([QgsPointXY(pt.x(), pt.y()) for pt in part])
				length = part_geom.length()
				if length <= 0.0:
					continue

				# Prepare chain point lists per part
				# Continuous chains per side: edges (left/right) and C points (left/right)
				chain_edge_left = []
				chain_edge_right = []
				chain_c_left = []
				chain_c_right = []

				# Stations: multiples of step_m; plotting requires exact user-specified matches only
				# We ensure the last station coincides with the part end so edges are emitted there too.
				stations = []
				if step_m > 0:
					n_steps = int(math.floor(length / step_m + 1e-9))
					for k in range(n_steps + 1):
						stations.append(round(k * step_m, 6))
				else:
					stations = [0.0]
				if not stations or abs(stations[-1] - length) > 1e-6:
					stations.append(float(round(length, 6)))
				stations = sorted(set(stations))

				# Mid-station for this part (nearest discrete station)
				# Used to pick a default SVG station when the user did not request specific ones.
				mid_station = None
				if stations:
					mid_station = min(stations, key=lambda s: abs(s - (length / 2.0)))

				# Small delta for tangent
				# We use a short segment around the station to compute a robust tangent direction.
				eps = min(0.5, step_m / 2.0) if step_m > 0 else 0.25

				for station in stations:
					p_center_pt = interp_point(part_geom, station)
					if p_center_pt is None:
						continue
					p_center = p_center_pt

					p_prev_pt = interp_point(part_geom, max(0.0, station - eps))
					p_next_pt = interp_point(part_geom, min(length, station + eps))
					if p_prev_pt is None or p_next_pt is None:
						continue
					tx = p_next_pt.x() - p_prev_pt.x()
					ty = p_next_pt.y() - p_prev_pt.y()
					norm = math.hypot(tx, ty)
					if norm == 0.0:
						continue
					tx /= norm
					ty /= norm
					nx = ty
					ny = -tx

					# Keep normal fixed to line direction: +s points to geometric right of the line
					# (No flipping by terrain; we decide uphill later from DEM at edges.)

					# Axis elevation z0 from DEM with global DEM offset, plus dz_axis
					# z0 is the trail surface level at the centerline.
					z0 = sample_dem(p_center)
					if z0 is not None:
						z0 = float(z0) + float(dem_offset) + float(dz_axis)
					if z0 is None:
						# Emit minimal points even if DEM sampling fails at axis (debug/visibility)
						ptA = QgsFeature(pt_fields)
						ptA.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(p_center.x(), p_center.y())))
						# Axis point: populate z fields with z0, keep area None
						ptA.setAttributes([station, "axis_center", "center", "edge", None, float(round(z0, 3)), float(round(z0, 3)), float(round(0.0, 3)), 0, area_cut_total, area_fill_total])
						if points_sink is not None and points_sink.addFeature(ptA, QgsFeatureSink.FastInsert):
							points_added += 1
						# trail edges in s for width-only placement
						s_L_dbg = -W / 2.0
						s_R_dbg = +W / 2.0
						for side, s_edge in [("left", s_L_dbg), ("right", s_R_dbg)]:
							xe = p_center.x() + nx * s_edge
							ye = p_center.y() + ny * s_edge
							pt = QgsFeature(pt_fields)
							ptype = "edge_left" if side == "left" else "edge_right"
							pt.setAttributes([station, ptype, side, "edge", None, None, None, None, 0])
							if points_sink is not None and points_sink.addFeature(pt, QgsFeatureSink.FastInsert):
								points_added += 1
						# Skip intersections/areas for this station
						continue

					# Build base DEM profile at DEM pixel spacing (piecewise-linear base)
					s_base = []
					z_base = []           # offset DEM (used everywhere)
					z_base_orig = []      # original DEM (for SVG reference only)
					s_cur = -R
					while s_cur <= R + 1e-9:
						x = p_center.x() + nx * s_cur
						y = p_center.y() + ny * s_cur
						z = sample_dem(QgsPointXY(x, y))
						s_base.append(float(round(s_cur, 3)))
						z_base_orig.append(None if z is None else float(round(z, 3)))
						# Apply global DEM offset for computations
						z_base.append(None if z is None else float(round(float(z) + float(dem_offset), 3)))
						s_cur += dem_step

					# Do not emit profile line features
					# We only keep arrays for computations and the optional plot.

					# (Axis center point will be emitted later with area attributes.)

					# trail edges in s (with +s pointing to right of line)
					# Left is negative, right is positive by our convention.
					s_L = -W / 2.0  # left edge
					s_R = +W / 2.0  # right edge
					# Helper to get DEM z at arbitrary s via linear interp on sampled profile (offset DEM)
					def z_dem_at_s(s_query):
						for i in range(len(s_base) - 1):
							s0, s1 = s_base[i], s_base[i + 1]
							if (s0 <= s_query <= s1) or (s1 <= s_query <= s0):
								z0v = z_base[i]
								z1v = z_base[i + 1]
								if z0v is None or z1v is None:
									return None
								if s1 == s0:
									return float(z0v)
								mu = (s_query - s0) / (s1 - s0)
								return float(round(z0v + mu * (z1v - z0v), 3))
						return None

					# Original DEM interpolation (for SVG reference)
					def z_dem_orig_at_s(s_query):
						for i in range(len(s_base) - 1):
							s0, s1 = s_base[i], s_base[i + 1]
							if (s0 <= s_query <= s1) or (s1 <= s_query <= s0):
								z0v = z_base_orig[i]
								z1v = z_base_orig[i + 1]
								if z0v is None or z1v is None:
									return None
								if s1 == s0:
									return float(z0v)
								mu = (s_query - s0) / (s1 - s0)
								return float(round(z0v + mu * (z1v - z0v), 3))
						return None

					zL_dem = z_dem_at_s(s_L)
					zR_dem = z_dem_at_s(s_R)

					# Optional override of axis elevation: full cut/fill sets z0 to edge min/max
					try:
						if axis_full_cut or axis_full_fill:
							cand = []
							if zL_dem is not None:
								cand.append(float(zL_dem))
							if zR_dem is not None:
								cand.append(float(zR_dem))
							if cand:
								z0 = (min(cand) if axis_full_cut else max(cand))
					except Exception:
						pass

					# Re-anchor z0 to the profile-sampled DEM at s=0 to ensure dz applies consistently
					# Only when not using full cut/fill axis override
					if not (axis_full_cut or axis_full_fill):
						z_axis_profile = z_dem_at_s(0.0)
						if z_axis_profile is not None:
							# z_axis_profile already includes dem_offset; add dz_axis
							z0 = float(z_axis_profile) + float(dz_axis)

					# Intersection finder that tests both up and down slopes from each edge.
					tanb = math.tan(beta1_rad)
					def _outward_distance(sq, edge_s):
						return (edge_s - sq) if sq <= edge_s else (sq - edge_s)

					def _bracket_index_for(edge_s):
						for i in range(len(s_base) - 1):
							if (s_base[i] <= edge_s <= s_base[i + 1]) or (s_base[i + 1] <= edge_s <= s_base[i]):
								return i
						return max(0, min(len(s_base) - 2, 0))

					def find_first_intersection(edge_s, side):
						idx0 = _bracket_index_for(edge_s)
						# Helper: test a single segment [sa,sb] with DEM end z's, return nearest hit if any
						def _test_segment(sa, sb, z_sa, z_sb):
							def z_up(sq):
								return z0 + tanb * _outward_distance(sq, edge_s)
							def z_dn(sq):
								return z0 - tanb * _outward_distance(sq, edge_s)
							f0_up = z_sa - z_up(sa)
							f1_up = z_sb - z_up(sb)
							f0_dn = z_sa - z_dn(sa)
							f1_dn = z_sb - z_dn(sb)
							cands = []
							def _add(mu, mode):
								ss = sa + mu * (sb - sa)
								zz = z_sa + mu * (z_sb - z_sa)
								cands.append((_outward_distance(ss, edge_s), float(round(ss, 3)), float(round(zz, 3)), mode))
							if f0_up == 0.0:
								_add(0.0, "up")
							if f1_up == 0.0:
								_add(1.0, "up")
							if f0_dn == 0.0:
								_add(0.0, "down")
							if f1_dn == 0.0:
								_add(1.0, "down")
							if f0_up * f1_up < 0.0:
								mu = -f0_up / (f1_up - f0_up)
								_add(mu, "up")
							if f0_dn * f1_dn < 0.0:
								mu = -f0_dn / (f1_dn - f0_dn)
								_add(mu, "down")
							return min(cands, key=lambda t: t[0]) if cands else None

						# 1) First partial segment from edge to nearest base breakpoint on outward side
						if side == "left":
							# outward decreasing s: neighbor is s_base[idx0]
							sa0 = edge_s
							sb0 = s_base[idx0]
							if sb0 < sa0:
								z_sa0 = z_dem_at_s(sa0)
								z_sb0 = z_base[idx0] if z_base[idx0] is not None else z_dem_at_s(sb0)
								best = _test_segment(sa0, sb0, z_sa0, z_sb0) if (z_sa0 is not None and z_sb0 is not None) else None
								if best:
									return best[1], best[2], best[3]
							# march further left over full base segments
							for i in range(idx0 - 1, -1, -1):
								sa = s_base[i + 1]
								sb = s_base[i]
								z_sa = z_base[i + 1] if z_base[i + 1] is not None else z_dem_at_s(sa)
								z_sb = z_base[i] if z_base[i] is not None else z_dem_at_s(sb)
								best = _test_segment(sa, sb, z_sa, z_sb) if (z_sa is not None and z_sb is not None) else None
								if best:
									return best[1], best[2], best[3]
						else:
							# right side: outward increasing s, neighbor is s_base[idx0+1]
							sa0 = edge_s
							sb0 = s_base[idx0 + 1] if (idx0 + 1) < len(s_base) else None
							if sb0 is not None and sb0 > sa0:
								z_sa0 = z_dem_at_s(sa0)
								z_sb0 = z_base[idx0 + 1] if z_base[idx0 + 1] is not None else z_dem_at_s(sb0)
								best = _test_segment(sa0, sb0, z_sa0, z_sb0) if (z_sa0 is not None and z_sb0 is not None) else None
								if best:
									return best[1], best[2], best[3]
							# march further right over full base segments
							for i in range(idx0 + 1, len(s_base) - 1):
								sa = s_base[i]
								sb = s_base[i + 1]
								z_sa = z_base[i] if z_base[i] is not None else z_dem_at_s(sa)
								z_sb = z_base[i + 1] if z_base[i + 1] is not None else z_dem_at_s(sb)
								best = _test_segment(sa, sb, z_sa, z_sb) if (z_sa is not None and z_sb is not None) else None
								if best:
									return best[1], best[2], best[3]

						return None, None, None

					# Compute intersections and emit points/segments
					sC_l, zC_l, dir_left = find_first_intersection(s_L, "left")
					sC_r, zC_r, dir_right = find_first_intersection(s_R, "right")

					# Design elevation function using chosen directions for each side
					def z_design_at_s(sq):
						s_lo = s_L if s_L <= s_R else s_R
						s_hi = s_R if s_R >= s_L else s_L
						if s_lo <= sq <= s_hi:
							return z0
						if sq <= s_L:
							d = (s_L - sq)
							return z0 + (tanb * d if dir_left == "up" else -tanb * d)
						if sq >= s_R:
							d = (sq - s_R)
							return z0 + (tanb * d if dir_right == "up" else -tanb * d)
						return z0

					c_points = []  # collect C points to emit after areas are known
					def emit_C(side, sC, zC, s_edge, slope_dir):
						nonlocal chain_c_left, chain_c_right, points_added
						# Create outputs tied to a found intersection C on one side:
						# - a short segment from the trail edge to C
						# - update simple side chains (left/right)
						if sC is None or zC is None:
							return
						xc = p_center.x() + nx * sC
						yc = p_center.y() + ny * sC
						# Collect C data for later point emission with area
						c_points.append({
							"side": side,
							"sC": sC,
							"zC": zC,
							"s_edge": s_edge,
							"slope_dir": slope_dir,
							"xc": xc,
							"yc": yc,
							"slope_run": abs(sC - s_edge),
						})

						# Add to C chains by side only
						if side == "left":
							chain_c_left.append(QgsPointXY(xc, yc))
						else:
							chain_c_right.append(QgsPointXY(xc, yc))

						# Segment from edge to C
						xe = p_center.x() + nx * s_edge
						ye = p_center.y() + ny * s_edge
						seg = QgsFeature(seg_fields)
						seg.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(xe, ye), QgsPointXY(xc, yc)]))
						dz_bc = float(round(zC - z0, 3))
						seg.setAttributes([station, side, "cut" if slope_dir == "up" else "fill", float(round(abs(sC - s_edge), 3)), dz_bc, "bc_line", None, _plot_name])
						segments_sink.addFeature(seg, QgsFeatureSink.FastInsert)

					emit_C("right", sC_r, zC_r, s_R, dir_right)
					emit_C("left", sC_l, zC_l, s_L, dir_left)

					# If a C point is missing on a side, add the edge point into the C chain to avoid gaps
					if sC_l is None:
						xe_fallback = p_center.x() + nx * s_L
						ye_fallback = p_center.y() + ny * s_L
						chain_c_left.append(QgsPointXY(xe_fallback, ye_fallback))
					if sC_r is None:
						xe_fallback = p_center.x() + nx * s_R
						ye_fallback = p_center.y() + ny * s_R
						chain_c_right.append(QgsPointXY(xe_fallback, ye_fallback))

					# Areas per station (computed for this station)
					# We integrate the vertical gap (DEM vs design) along s using the trapezoid rule
					# to get cross-sectional areas of cut and fill. Units are square meters per station.
					def build_samples(s_start, s_end):
						if s_start is None or s_end is None:
							return []
						if s_start > s_end:
							s_start, s_end = s_end, s_start
						ss = []
						ss.append(s_start)
						for sv in s_base:
							if s_start < sv < s_end:
								ss.append(sv)
						ss.append(s_end)
						ss = sorted(set(round(x, 6) for x in ss))
						out = []
						for sqq in ss:
							zdem = z_dem_at_s(sqq)
							if zdem is None:
								continue
							zdes = z_design_at_s(sqq)
							out.append((float(sqq), float(zdem), float(zdes)))
						return out

					def trapz_area(samples, mode):
						# Trapezoidal integration of positive gaps only.
						# mode='cut' uses (DEM - design) if positive, 'fill' uses (design - DEM).
						if not samples or len(samples) < 2:
							return None
						acc = 0.0
						for i in range(len(samples) - 1):
							s0, z0d, z0s = samples[i]
							s1, z1d, z1s = samples[i + 1]
							g0 = (z0d - z0s) if mode == 'cut' else (z0s - z0d)
							g1 = (z1d - z1s) if mode == 'cut' else (z1s - z1d)
							g0 = max(0.0, g0)
							g1 = max(0.0, g1)
							ds = (s1 - s0)
							acc += 0.5 * (g0 + g1) * ds
						return float(round(acc, 3))

					# Trail area between edges
					trail_samples = build_samples(s_L, s_R)
					area_cut_trail = trapz_area(trail_samples, 'cut') if trail_samples else None
					area_fill_trail = trapz_area(trail_samples, 'fill') if trail_samples else None

					# Areas on each side: edge to C => total cut/fill only
					u_samples = build_samples(s_L, sC_l) if sC_l is not None else []
					d_samples = build_samples(s_R, sC_r) if sC_r is not None else []
					area_cut_up = trapz_area(u_samples, 'cut') if u_samples else None
					area_fill_up = trapz_area(u_samples, 'fill') if u_samples else None
					area_cut_down = trapz_area(d_samples, 'cut') if d_samples else None
					area_fill_down = trapz_area(d_samples, 'fill') if d_samples else None

					cut_vals = [v for v in [area_cut_trail, area_cut_up, area_cut_down] if v is not None]
					fill_vals = [v for v in [area_fill_trail, area_fill_up, area_fill_down] if v is not None]
					area_cut_total = float(round(sum(cut_vals), 3)) if cut_vals else None
					area_fill_total = float(round(sum(fill_vals), 3)) if fill_vals else None

					# Assign signed per-point areas to output features for the four intervals
					# Convention: cut positive, fill negative; axis receives None
					def _signed_area(cut_val, fill_val):
						c = cut_val if isinstance(cut_val, (int, float)) else 0.0
						f = fill_val if isinstance(fill_val, (int, float)) else 0.0
						return float(round(c - f, 3)) if (c or f) else None

					# Per-side trail areas from axis to each edge (signed)
					left_trail_samples = build_samples(0.0, s_L)
					right_trail_samples = build_samples(0.0, s_R)
					left_trail_area = _signed_area(
						trapz_area(left_trail_samples, 'cut') if left_trail_samples else None,
						trapz_area(left_trail_samples, 'fill') if left_trail_samples else None
					)
					right_trail_area = _signed_area(
						trapz_area(right_trail_samples, 'cut') if right_trail_samples else None,
						trapz_area(right_trail_samples, 'fill') if right_trail_samples else None
					)

					# Left slope: edge to C_left
					left_slope_area = _signed_area(area_cut_up, area_fill_up)
					# Right slope: edge to C_right
					right_slope_area = _signed_area(area_cut_down, area_fill_down)

					# Determine if this station is selected for SVG (for svg_sel flag on points)
					svg_flag = 1 if (svg_all or (requested_stations and any(abs(station - req) <= 1e-6 for req in requested_stations))) else 0

					# Emit points now with proper area fields
					# Axis center point
					ptA = QgsFeature(pt_fields)
					ptA.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(p_center.x(), p_center.y())))
					ptA.setAttributes([station, "axis_center", "center", "edge", None, float(round(z0, 3)), float(round(z0, 3)), float(round(0.0, 3)), svg_flag, area_cut_total, area_fill_total])
					if points_sink is not None and points_sink.addFeature(ptA, QgsFeatureSink.FastInsert):
						points_added += 1
					# Edge points with per-side area from axis to that edge (signed)
					xeL, yeL = p_center.x() + nx * s_L, p_center.y() + ny * s_L
					xeR, yeR = p_center.x() + nx * s_R, p_center.y() + ny * s_R
					# Append to edge chains
					chain_edge_left.append(QgsPointXY(xeL, yeL))
					chain_edge_right.append(QgsPointXY(xeR, yeR))
					for side, s_edge, xe, ye, edge_dem, area_part in [
						("left", s_L, xeL, yeL, zL_dem, left_trail_area),
						("right", s_R, xeR, yeR, zR_dem, right_trail_area),
					]:
						pt = QgsFeature(pt_fields)
						pt.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(xe, ye)))
						ptype = "edge_left" if side == "left" else "edge_right"
						z_val = float(round(z0, 3))
						edge_kind = "edge"
						if edge_dem is not None:
							edge_kind = "cut" if float(edge_dem) > float(z0) else ("fill" if float(edge_dem) < float(z0) else "edge")
						# Assign areas to separate fields: area_cut (positive), area_fill (negative)
						area_cut_val = area_part if (isinstance(area_part, (int, float)) and area_part is not None and area_part > 0) else None
						area_fill_val = area_part if (isinstance(area_part, (int, float)) and area_part is not None and area_part < 0) else None
						pt.setAttributes([station, ptype, side, edge_kind, None, z_val, float(round(z0, 3)), float(round(z_val - z0, 3)), svg_flag, area_cut_total, area_fill_total])
						if points_sink is not None and points_sink.addFeature(pt, QgsFeatureSink.FastInsert):
							points_added += 1
					# C points with slope areas
					for cp in c_points:
						area_part = left_slope_area if cp["side"] == "left" else right_slope_area
						if area_part is None:
							continue
						ptype = "slope_right" if cp["side"] == "right" else "slope_left"
						# Derive cut/fill from area sign for consistency
						tri_kind = "cut" if (isinstance(area_part, (int, float)) and area_part >= 0.0) else "fill"
						z_val = float(round(cp["zC"], 3))
						ptC = QgsFeature(pt_fields)
						ptC.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(cp["xc"], cp["yc"])) )
						area_cut_val = area_part if (isinstance(area_part, (int, float)) and area_part is not None and area_part > 0) else None
						area_fill_val = area_part if (isinstance(area_part, (int, float)) and area_part is not None and area_part < 0) else None
						ptC.setAttributes([station, ptype, cp["side"], tri_kind, float(round(cp["slope_run"], 3)), z_val, float(round(z0, 3)), float(round(z_val - z0, 3)), svg_flag, area_cut_total, area_fill_total])
						if points_sink is not None and points_sink.addFeature(ptC, QgsFeatureSink.FastInsert):
							points_added += 1

					# SVG preview logic per-station
					# If enabled, we draw an equal-scale cross-section with key markers and labels.
					should_plot = False
					filename_suffix = None
					if True:
						# Tolerance for matching requested stations
						tol = max(step_m / 2.0, 1e-6)
						if svg_all and svg_count < MAX_SVGS:
							should_plot = True
							filename_suffix = f"s_{int(round(station))}"
						elif requested_stations and svg_count < MAX_SVGS:
							for req in requested_stations:
								if abs(station - req) <= 1e-6:
									should_plot = True
									filename_suffix = f"s_{int(round(req))}"
									matched_stations_for_feature.add(req)
									break
						else:
							# If no requested stations and not all, or cap reached, skip SVG entirely
							should_plot = False

					if should_plot:
						try:
							# Debug: log sign conventions for the plotted station
							sCl_str = "None" if sC_l is None else f"{sC_l:.3f}"
							sCr_str = "None" if sC_r is None else f"{sC_r:.3f}"
							feedback.pushInfo(
								f"SVG DBG station={station:.2f}: s_L={s_L:.3f}, s_R={s_R:.3f}, sC_l={sCl_str}, sC_r={sCr_str}, dir_left={dir_left}, dir_right={dir_right}"
							)
							# Build arrays for plotting
							# Reuse DEM base breakpoints for plotting (exact piecewise-linear DEM)
							plot_s = [float(sv) for sv in s_base]
							plot_dem = [None if zv is None else float(zv) for zv in z_base]

							# Create SVG
							# We show: DEM profile (smoothed for gaps), flat trail between edges,
							# side slopes to C points, helper dimensions, and labels A/L/R/C.
							import matplotlib
							matplotlib.use("Agg")
							from matplotlib import pyplot as plt
							from matplotlib.ticker import MultipleLocator

							fig, ax = plt.subplots(1, 1, figsize=(8, 5))
							# Filter None by simple forward-fill for plotting continuity
							dem_vals = []
							last_ok = None
							for v in plot_dem:
								if v is None:
									dem_vals.append(last_ok)
								else:
									dem_vals.append(v)
									last_ok = v

							# Trail surface between edges (left to right)
							trail_line = ax.plot([s_L, s_R], [z0, z0], label="Planum", color="#000000", lw=1.0)[0]
							# Single continuous DEM line (offset DEM used for computations)
							if abs(dem_offset) <= 1e-9:
								_dem_label = "DEM"
								_dem_color = "#2ca02c"
							else:
								_dem_label = f"DEM (dz={dem_offset:.2f} m)"
								_dem_color = "#8B4513"
							dem_line = ax.plot(plot_s, dem_vals, color=_dem_color, ls="--" , lw=1.0, label=_dem_label)[0]

							# Road/trail surface (trapezoid above Planum)
							# Compute isosceles trapezoid with inward-tilted sides; horizontal base/top
							if surf_height_eff > 0.0:
								try:
									# Horizontal base: (s_L,z0) -> (s_R,z0)
									# Horizontal top at z0+H; inward tilt angle from vertical defines horizontal inset per side
									H = float(surf_height_eff)
									# Angle is from horizontal; horizontal inset per side = H / tan(angle)
									if vertical_sides:
										inset_per_side = 0.0
									else:
										_tanv = math.tan(surf_xfall_rad)
										inset_per_side = (H / _tanv) if _tanv != 0.0 else 0.0
									# Clamp inset so top does not invert
									max_inset = max(0.0, (s_R - s_L) / 2.0 - 1e-6)
									inset = min(max_inset, max(0.0, inset_per_side))
									top_z = z0 + H
									top_L = s_L + inset
									top_R = s_R - inset
									# Determine top crossfall per user mode (visual only)
									# 0: Flat; 1: Inward (toward uphill, +tilt°); 2: Outward (toward downhill, +tilt°)
									g_top = 0.0
									try:
										zL_dem = z_dem_at_s(s_L)
										zR_dem = z_dem_at_s(s_R)
										dzLR = None if (zL_dem is None or zR_dem is None) else (float(zR_dem) - float(zL_dem))
										if surf_top_mode == 1:  # inward: tilt to uphill side
											# uphill is the higher DEM edge; slope down toward lower edge
											if dzLR is not None:
												# Positive dzLR => right is uphill; negative => left is uphill; zero => right
												dir_sign = -1.0 if dzLR > 0 else (1.0 if dzLR < 0 else -1.0)
												g_top = math.tan(math.radians(max(0, min(7, surf_top_tilt_deg)))) * dir_sign
										elif surf_top_mode == 2:  # outward: tilt to downhill side
											if dzLR is not None:
												dir_sign = 1.0 if dzLR > 0 else (-1.0 if dzLR < 0 else 1.0)
												g_top = math.tan(math.radians(max(0, min(7, surf_top_tilt_deg)))) * dir_sign
									except Exception:
										g_top = 0.0

									# Build top edge y with crossfall, rotating around base midpoint so average height = top_z
									s_mid = 0.5 * (s_L + s_R)
									# Offsets from pivot
									dL = (top_L - s_mid)
									dR = (top_R - s_mid)
									# To preserve constant height H, symmetrically offset around top_z
									top_L_y = top_z + g_top * dL
									top_R_y = top_z + g_top * dR
									# Draw top outline (no legend)
									ax.plot([top_L, top_R], [top_L_y, top_R_y], color="#666666", lw=1.0, label="_nolegend_", zorder=4)
									# Label tilt angle q at horizontal center, parallel to surface top
									try:
										if surf_top_mode in (1, 2) and surf_top_tilt_deg > 0:
											theta = math.atan2((top_R_y - top_L_y), max(1e-9, (top_R - top_L)))
											xc = 0.5 * (top_L + top_R)
											yc = 0.5 * (top_L_y + top_R_y)
											# offset slightly along surface normal so text sits above the top edge
											nx_txt = -math.sin(theta)
											ny_txt =  math.cos(theta)
											offset_d = max(0.05, H * 0.005)
											x_txt = xc + nx_txt * offset_d
											y_txt = yc + ny_txt * offset_d
											ax.text(x_txt, y_txt, f"- q: {int(max(0, min(7, surf_top_tilt_deg)))}° -", rotation=math.degrees(theta), rotation_mode="anchor", ha="center", va="bottom", fontsize=8, color="#666666", bbox=MPL_BBOX, zorder=1000)
									except Exception:
										pass
									# Filled polygon: base left -> base right -> top right -> top left (with tilt)
									# Compute area and label text for legend (shoelace for the quadrilateral)
									try:
										# vertices in order: base left -> base right -> top right -> top left
										xs = [s_L, s_R, top_R, top_L]
										ys = [z0, z0, top_R_y, top_L_y]
										acc = 0.0
										for i in range(4):
											j = (i + 1) % 4
											acc += xs[i] * ys[j] - xs[j] * ys[i]
										area_top = float(round(abs(acc) * 0.5, 2))
									except Exception:
										# fallback to trapezoid area if any issue
										width_base = (s_R - s_L)
										width_top = max(0.0, (top_R - top_L))
										area_top = float(round(H * (width_base + width_top) / 2.0, 2))
									surface_label_text = f"Fahrbahnkörper ($h_{{A}}$={float(round(H, 2)):.2f} m, $A$={area_top:.2f} m²)"
									surface_poly = ax.fill([s_L, s_R, top_R, top_L], [z0, z0, top_R_y, top_L_y], color="#999999", alpha=0.4, label="_nolegend_", zorder=3)[0]
								except Exception:
									pass

							# Original DEM for reference and shaded band only if non-zero offset
							if abs(dem_offset) > 1e-9:
								try:
									import numpy as np
									x_np_ref = np.array([float(sv) for sv in plot_s], dtype=float)
									dem_np_off = np.array([z_dem_at_s(sv) for sv in x_np_ref], dtype=float)
									dem_np_orig = np.array([z_dem_orig_at_s(sv) for sv in x_np_ref], dtype=float)
									valid_ref = np.isfinite(dem_np_off) & np.isfinite(dem_np_orig)
									if np.any(valid_ref):
										orig_dem_line = ax.plot(x_np_ref[valid_ref], dem_np_orig[valid_ref], color="#2ca02c", ls=":", lw=1.0, label="DEM (Ursprung)")[0]
										dem_offset_band = ax.fill_between(x_np_ref[valid_ref], dem_np_orig[valid_ref], dem_np_off[valid_ref], color="#c7a27c", alpha=0.8, label="Oberboden")
								except Exception:
									pass

							# Shaded cut/fill areas between DEM and design surface
							try:
								import numpy as np
								# Start from DEM base breakpoints (exact PL fit)
								x_np = np.array([float(sv) for sv in plot_s], dtype=float)
								dem_np = np.array([z_dem_at_s(sv) for sv in x_np], dtype=float)
								des_np = np.array([z_design_at_s(sv) for sv in x_np], dtype=float)
								# Ensure critical breakpoints exist: 0, s_L, s_R, sC_l, sC_r
								critical = [0.0]
								if s_L is not None and (-R - 1e-9) <= s_L <= (R + 1e-9):
									critical.append(float(s_L))
								if s_R is not None and (-R - 1e-9) <= s_R <= (R + 1e-9):
									critical.append(float(s_R))
								if sC_l is not None and (-R - 1e-9) <= sC_l <= (R + 1e-9):
									critical.append(float(sC_l))
								if sC_r is not None and (-R - 1e-9) <= sC_r <= (R + 1e-9):
									critical.append(float(sC_r))
								for sx in critical:
									if not np.any(np.isclose(x_np, sx, atol=1e-9)):
										x_np = np.append(x_np, sx)
										dem_val = z_dem_at_s(sx)
										des_val = z_design_at_s(sx)
										dem_np = np.append(dem_np, (np.nan if dem_val is None else float(dem_val)))
										des_np = np.append(des_np, (np.nan if des_val is None else float(des_val)))
								# Sort by s
								ord_idx = np.argsort(x_np)
								x_np = x_np[ord_idx]
								dem_np = dem_np[ord_idx]
								des_np = des_np[ord_idx]
								valid = np.isfinite(dem_np) & np.isfinite(des_np)
								# Masks limited to axis→C per side
								left_mask = None
								right_mask = None
								if sC_l is not None:
									lo = min(0.0, sC_l)
									hi = max(0.0, sC_l)
									left_mask = (x_np >= lo) & (x_np <= hi)
								if sC_r is not None:
									lo = min(0.0, sC_r)
									hi = max(0.0, sC_r)
									right_mask = (x_np >= lo) & (x_np <= hi)
								is_cut = valid & (dem_np > des_np)
								is_fill = valid & (des_np > dem_np)
								# Compose masks
								mask_cut = None
								mask_fill = None
								if left_mask is not None:
									mask_cut = left_mask & is_cut if mask_cut is None else (mask_cut | (left_mask & is_cut))
									mask_fill = left_mask & is_fill if mask_fill is None else (mask_fill | (left_mask & is_fill))
								if right_mask is not None:
									mask_cut = right_mask & is_cut if mask_cut is None else (mask_cut | (right_mask & is_cut))
									mask_fill = right_mask & is_fill if mask_fill is None else (mask_fill | (right_mask & is_fill))
								label_cut = f"Einschnitt ($A$={area_cut_total:.2f} m²)" if area_cut_total is not None else "Einschnitt"
								label_fill = f"Auftrag ($A$={area_fill_total:.2f} m²)" if area_fill_total is not None else "Auftrag"
								cut_poly = fill_poly = None
								if mask_cut is not None and np.any(mask_cut):
									cut_poly = ax.fill_between(x_np, dem_np, des_np, where=mask_cut, interpolate=True, color="#f4cccc", alpha=0.6, label=label_cut, zorder=1)
								if mask_fill is not None and np.any(mask_fill):
									fill_poly = ax.fill_between(x_np, dem_np, des_np, where=mask_fill, interpolate=True, color="#d9ead3", alpha=0.6, label=label_fill, zorder=1)
							except Exception:
								pass

							# Slope lines from edges to C (if present)
							_added_slope_legend = False
							slope_handle = None
							if sC_l is not None and zC_l is not None:
								_h = ax.plot([s_L, sC_l], [z0, zC_l], label=("Böschung" if not _added_slope_legend else "_nolegend_"), color="#cc6677")[0]
								if not _added_slope_legend:
									slope_handle = _h
								_added_slope_legend = True
							if sC_r is not None and zC_r is not None:
								_h = ax.plot([s_R, sC_r], [z0, zC_r], label=("Böschung" if not _added_slope_legend else "_nolegend_"), color="#cc6677")[0]
								if not _added_slope_legend and slope_handle is None:
									slope_handle = _h
								_added_slope_legend = True

							# Helper lines and labels
							if sC_l is not None and zC_l is not None:
								ax.plot([0.0, sC_l], [z0, z0], color="#808080", ls=":", lw=1.0, zorder=-1)
								ax.plot([sC_l, sC_l], [z0, zC_l], color="#808080", ls=":", lw=1.0, zorder=1.1)
								dzL = (zC_l - z0)
								y_mid_L = z0 + (zC_l - z0) / 2.0
								_dx = 0.08
								if sC_l < 0:
									ax.text(sC_l - _dx, y_mid_L, f"h: {dzL:.2f} m", color="#808080", fontsize=8, va="center", ha="right", bbox=MPL_BBOX, zorder=1000)
								else:
									ax.text(sC_l + _dx, y_mid_L, f"h: {dzL:.2f} m", color="#808080", fontsize=8, va="center", ha="left", bbox=MPL_BBOX, zorder=1000)
								# Base run label (right of C-left)
								run_L = abs(sC_l - s_L)
								ymin_axis_lbl = ax.get_ylim()[0]-0.5
								ax.text(sC_l + 0.08, ymin_axis_lbl, f"w: {run_L:.2f} m", color="#808080", fontsize=8, va="bottom", ha="left", bbox=MPL_BBOX, zorder=1000)

							if sC_r is not None and zC_r is not None:
								ax.plot([0.0, sC_r], [z0, z0], color="#808080", ls=":", lw=1, zorder=-1)
								ax.plot([sC_r, sC_r], [z0, zC_r], color="#808080", ls=":", lw=1, zorder=1.1)
								dzR = (zC_r - z0)
								y_mid_R = z0 + (zC_r - z0) / 2.0
								_dx = 0.08
								if sC_r < 0:
									ax.text(sC_r - _dx, y_mid_R, f"h: {dzR:.2f} m", color="#808080", fontsize=8, va="center", ha="right", bbox=MPL_BBOX, zorder=1000)
								else:
									ax.text(sC_r + _dx, y_mid_R, f"h: {dzR:.2f} m", color="#808080", fontsize=8, va="center", ha="left", bbox=MPL_BBOX, zorder=1000)
								# Base run label (left of C-right)
								run_R = abs(sC_r - s_R)
								ymin_axis_lbl = ax.get_ylim()[0]-0.5
								ax.text(sC_r - 0.08, ymin_axis_lbl, f"w: {run_R:.2f} m", color="#808080", fontsize=8, va="bottom", ha="right", bbox=MPL_BBOX, zorder=1000)

							# Markers (reduced sizes)
							ax.scatter([0.0], [z0], color="#000000", s=7, zorder=5)
							ax.text(0.0, z0, " A", fontsize=8, va="bottom", zorder=1000)
							# Edge colors by cut/fill vs axis z0
							edge_color_L = ("#cc3333" if float(zL_dem) > float(z0) else ("#2ca02c" if float(zL_dem) < float(z0) else "#000000"))
							if zR_dem is not None:
								edge_color_R = ("#cc3333" if float(zR_dem) > float(z0) else ("#2ca02c" if float(zR_dem) < float(z0) else "#000000"))
							ax.scatter([s_L], [z0], color=edge_color_L, s=6, zorder=5)
							ax.scatter([s_R], [z0], color=edge_color_R, s=6, zorder=5)
							ax.text(s_L, z0, " L", fontsize=8, va="bottom", zorder=1000)
							ax.text(s_R, z0, " R", fontsize=8, va="bottom", zorder=1000)
							# Classify C markers by cut/fill using area sign first, fallback to z comparison
							if sC_l is not None and zC_l is not None:
								tri_left = None
								if 'left_slope_area' in locals() and left_slope_area is not None:
									tri_left = "cut" if float(left_slope_area) >= 0.0 else "fill"
								else:
									try:
										tri_left = "cut" if float(zC_l) > float(z_design_at_s(sC_l)) else "fill"
									except Exception:
										tri_left = None
								color_left_C = "#cc3333" if tri_left == "cut" else ("#2ca02c" if tri_left == "fill" else "#117733")
								ax.scatter([sC_l], [zC_l], color=color_left_C, s=8, zorder=6)
								# ax.text(sC_l, zC_l, " C_left", fontsize=8, va="bottom", color="#117733", bbox=dict(facecolor="white", edgecolor="none", pad=0.1, alpha=0.7))
							if sC_r is not None and zC_r is not None:
								tri_right = None
								if 'right_slope_area' in locals() and right_slope_area is not None:
									tri_right = "cut" if float(right_slope_area) >= 0.0 else "fill"
								else:
									try:
										tri_right = "cut" if float(zC_r) > float(z_design_at_s(sC_r)) else "fill"
									except Exception:
										tri_right = None
								color_right_C = "#cc3333" if tri_right == "cut" else ("#2ca02c" if tri_right == "fill" else "#117733")
								ax.scatter([sC_r], [zC_r], color=color_right_C, s=8, zorder=6)
								# ax.text(sC_r, zC_r, " C_right", fontsize=8, va="bottom", color="#117733", bbox=dict(facecolor="white", edgecolor="none", pad=0.1, alpha=0.7))
							# Highlight A with a white circle overlay
							ax.scatter([0.0], [z0], facecolors="none", edgecolors="#000000", s=20, linewidths=0.8, zorder=6)

							# If axis z was overridden (dz != 0 or full cut/fill), draw vertical dashed line from DEM to A with height label
							try:
								_z_dem_at_axis = z_dem_at_s(0.0)
								if _z_dem_at_axis is not None:
									_dhA = float(z0) - float(_z_dem_at_axis)
									if (abs(_dhA) > 1e-9) and (abs(dz_axis) > 1e-9 or axis_full_cut or axis_full_fill):
										_y0 = float(min(z0, _z_dem_at_axis))
										_y1 = float(max(z0, _z_dem_at_axis))
										ax.plot([0.0, 0.0], [_y0, _y1], color="#ff0000", ls="--", lw=1, zorder=-1)
										_ym = _y0 + (_y1 - _y0) / 2.0
										ax.text(0.08, _ym, f"h: {_dhA:.2f} m", color="#ff0000", fontsize=8, va="center", ha="left", bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.7), zorder=1000)
							except Exception:
								pass

							xlab = ax.set_xlabel("Querschnitt [m]")
							xlab.set_zorder(1000)
							ylab = ax.set_ylabel("Seehöhe [m]")
							ylab.set_zorder(1000)
							title_text = ax.set_title(f"Systemschnitt {_plot_name} — Station {station}")
							title_text.set_zorder(1000)
							ax.set_aspect("equal", adjustable="datalim")
							ax.xaxis.set_major_locator(MultipleLocator(1.0))
							ax.yaxis.set_major_locator(MultipleLocator(1.0))
							ax.grid(True, which="major", axis="both", alpha=0.3)
							# Make axes frame look like grid lines and shrink tick labels slightly
							try:
								import matplotlib as mpl
								grid_color = mpl.rcParams.get("grid.color", "#b0b0b0")
								grid_alpha = float(mpl.rcParams.get("grid.alpha", 0.3))
								grid_lw = float(mpl.rcParams.get("grid.linewidth", 0.8))
								for _sp in ax.spines.values():
									_sp.set_edgecolor(grid_color)
									_sp.set_alpha(grid_alpha)
									_sp.set_linewidth(grid_lw)
								# Reduce tick label font size a tad
								ax.tick_params(axis="both", labelsize=8)
							except Exception:
								pass
							# Ensure tick labels draw above other artists
							try:
								for _tick in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
									_tick.set_zorder(1000)
							except Exception:
								pass
							# Ordered legend: Böschung, cut, fill, DEM, Trail
							# DEM handle already has the correct label from plotting above
							handles = []
							labels = []
							if slope_handle is not None:
								handles.append(slope_handle)
								# Show Böschung with rise/run as reduced fraction from tan(beta_1)
								try:
									import fractions
									_tb = math.tan(beta1_rad)
									frac = fractions.Fraction(_tb).limit_denominator(10)
									labels.append(f"Böschung ({frac.numerator}:{frac.denominator})")
								except Exception:
									labels.append(f"Böschung ({beta1_deg:.1f}°)")
							if 'label_cut' in locals() and 'cut_poly' in locals() and cut_poly is not None:
								handles.append(cut_poly)
								labels.append(label_cut)
							if 'label_fill' in locals() and 'fill_poly' in locals() and fill_poly is not None:
								handles.append(fill_poly)
								labels.append(label_fill)
							handles.append(dem_line)
							labels.append(dem_line.get_label())
							# Include original DEM (surface) when offset is used
							if 'orig_dem_line' in locals():
								handles.append(orig_dem_line)
								labels.append(orig_dem_line.get_label())
							# Include the DEM offset band (Oberboden) when present
							if 'dem_offset_band' in locals():
								handles.append(dem_offset_band)
								labels.append(dem_offset_band.get_label())
							handles.append(trail_line)
							labels.append(trail_line.get_label())
							# Place Fahrbahnkörper last in legend if present
							if 'surface_poly' in locals() and surface_poly is not None:
								_lbl = surface_label_text if 'surface_label_text' in locals() else "Fahrbahnkörper"
								handles.append(surface_poly)
								labels.append(_lbl)
							leg = ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5))
							leg.set_zorder(1000)
							try:
								for t in leg.get_texts():
									lbl = t.get_text().lower().strip()
									if lbl.startswith("einschnitt"):
										t.set_color("#cc3333")
									elif lbl.startswith("auftrag"):
										t.set_color("#2ca02c")
							except Exception:
								pass

							# Vertical padding of 1 m on y-axis
							_y_vals = [v for v in dem_vals if v is not None]
							_y_vals.extend([z0])
							if sC_l is not None:
								_y_vals.append(zC_l)
							if sC_r is not None:
								_y_vals.append(zC_r)
							if _y_vals:
								ymin = min(_y_vals) - 1.0 # lower y-padding in plot area
								ymax = max(_y_vals) + 0.5 # upper y-padding in plot area
								ax.set_ylim(ymin, ymax)

								# Grey underlay from plot y-min to DEM (offset if active)
								try:
									import numpy as np
									y_base = ax.get_ylim()[0]
									x_under = np.array([float(sv) for sv in plot_s], dtype=float)
									y_under = np.array([z_dem_at_s(sv) for sv in x_under], dtype=float)
									finite_mask = np.isfinite(y_under)
									if np.any(finite_mask):
										ax.fill_between(x_under[finite_mask], y_base, y_under[finite_mask], color="#cccccc", alpha=0.35, label="_nolegend_", zorder=0)
								except Exception:
									pass

							# Enforce 1:1 axis scaling without forcing square area
							# Determine spans and set aspect equal (1:1 units)
							xmin, xmax = ax.get_xlim()
							ymin, ymax = ax.get_ylim()
							xspan = max(1e-9, (xmax - xmin))
							yspan = max(1e-9, (ymax - ymin))
							ax.set_aspect("equal", adjustable="datalim")
							# Compute figure size so that 1 m ~ svg_px_per_m pixels (independently for x and y)
							# Matplotlib fig size is in inches; 1 inch ~ 96 px in SVG backend
							px_per_in = 96.0
							fig_w_in = max(4.0, (xspan * svg_px_per_m) / px_per_in)
							fig_h_in = max(3.0, (yspan * svg_px_per_m) / px_per_in)
							fig.set_size_inches(fig_w_in, fig_h_in)

							# (Removed in-plot area annotation; areas shown in legend labels only)

							# Save to project directory (force overwrite via atomic replace)
							proj_file = QgsProject.instance().fileName() or ""
							proj_dir = os.path.dirname(proj_file) if proj_file else os.getcwd()
							if filename_suffix:
								out_name = f"section_{_plot_name}_{filename_suffix}.svg"
							else:
								out_name = f"section_{_plot_name}.svg"
							out_path = os.path.join(proj_dir, out_name)
							# Emit transect line feature for the plotted station: from s=-R to s=+R
							try:
								xL = p_center.x() + nx * (-R)
								yL = p_center.y() + ny * (-R)
								xR = p_center.x() + nx * (+R)
								yR = p_center.y() + ny * (+R)
								tran = QgsFeature(seg_fields)
								tran.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(xL, yL), QgsPointXY(xR, yR)]))
								tran.setAttributes([station, None, None, None, None, "transect_line", float(round(2.0 * R, 3)), _plot_name])
								segments_sink.addFeature(tran, QgsFeatureSink.FastInsert)
							except Exception:
								pass

							fig.tight_layout()
							# Test: dashed vertical lines at L and R clipped at planum (y from ymin to z0)
							ymin_axis = ax.get_ylim()[0]
							ax.vlines(s_L, ymin_axis, z0, colors="#808080", linestyles=":", linewidth=1.2, zorder=1.5, clip_on=False)
							ax.vlines(s_R, ymin_axis, z0, colors="#808080", linestyles=":", linewidth=1.2, zorder=1.5, clip_on=False)
							if sC_l is not None and zC_l is not None:
								ax.vlines(sC_l, ymin_axis, zC_l, colors="#808080", linestyles="-", linewidth=1.2, zorder=1.5, clip_on=False)
							if sC_r is not None and zC_r is not None:
								ax.vlines(sC_r, ymin_axis, zC_r, colors="#808080", linestyles="-", linewidth=1.2, zorder=1.5, clip_on=False)
							fig.tight_layout()
							tmp_path = out_path + ".tmp"
							try:
								fig.savefig(tmp_path, format="svg", bbox_inches="tight", pad_inches=0.1)
								os.replace(tmp_path, out_path)
								# count this SVG (open silently if requested)
								svg_count += 1
							finally:
								plt.close(fig)
							if svg_open:
								try:
									QDesktopServices.openUrl(QUrl.fromLocalFile(out_path))
								except Exception:
									pass
							svg_done_for_feature = True
						except Exception as _e:
							feedback.reportError(f"SVG generation failed: {_e}")

				# After stations of this part: emit simple chains if they have >=2 points
				if len(chain_edge_left) >= 2:
					f_el = QgsFeature(chain_fields2)
					f_el.setGeometry(QgsGeometry.fromPolylineXY(chain_edge_left))
					f_el.setAttributes(["left", "trail_edge"])
					chain_edges_sink.addFeature(f_el, QgsFeatureSink.FastInsert)
				if len(chain_edge_right) >= 2:
					f_er = QgsFeature(chain_fields2)
					f_er.setGeometry(QgsGeometry.fromPolylineXY(chain_edge_right))
					f_er.setAttributes(["right", "trail_edge"])
					chain_edges_sink.addFeature(f_er, QgsFeatureSink.FastInsert)
				# Emit C chains per side
				if len(chain_c_left) >= 2:
					f_cl = QgsFeature(chain_fields2)
					f_cl.setGeometry(QgsGeometry.fromPolylineXY(chain_c_left))
					f_cl.setAttributes(["left", "slope_edge"])
					chain_c_sink.addFeature(f_cl, QgsFeatureSink.FastInsert)
				if len(chain_c_right) >= 2:
					f_cr = QgsFeature(chain_fields2)
					f_cr.setGeometry(QgsGeometry.fromPolylineXY(chain_c_right))
					f_cr.setAttributes(["right", "slope_edge"])
					chain_c_sink.addFeature(f_cr, QgsFeatureSink.FastInsert)


			# If the user requested stations that were not matched, report it
			if requested_stations:
				for req in requested_stations:
					if req not in matched_stations_for_feature:
						feedback.pushInfo(f"SVG: Station {req} m not found on feature {feat.id()} (no plot).")

		# Basic styling: points size 2px, lines width 0.2, all white; highlight svg stations
		def _style_point_layer(layer_id):
			layer = context.temporaryLayerStore().mapLayer(layer_id) or QgsProject.instance().mapLayer(layer_id)
			if not layer:
				return
			# Pre-style by svg_sel and cut/fill: svg stations = white circle; cut=red; fill=green; edges/axis=black
			root = QgsRuleBasedRenderer.Rule(None)
			# svg station highlight (on top) — for non-axis points only
			rule_svg = QgsRuleBasedRenderer.Rule(QgsMarkerSymbol.createSimple({"color": "white", "outline_color": "black", "size": "1.1", "outline_width": "0.4", "name": "circle"}))
			rule_svg.setFilterExpression("\"svg_sel\" = 1 AND \"ptype\" <> 'axis_center'")
			rule_svg.setLabel("SVG-Stationen")
			root.appendChild(rule_svg)
			# cut/fill/edge rules with 50% size
			rule_cut = QgsRuleBasedRenderer.Rule(QgsMarkerSymbol.createSimple({"color": "#cc3333", "outline_color": "#cc3333", "size": "0.8", "outline_width": "0"}))
			rule_cut.setFilterExpression("\"ptype\" <> 'axis_center' AND \"cut_fill\" = 'cut'")
			rule_cut.setLabel("Einschnitt (cut)")
			root.appendChild(rule_cut)
			rule_fill = QgsRuleBasedRenderer.Rule(QgsMarkerSymbol.createSimple({"color": "#2ca02c", "outline_color": "#2ca02c", "size": "0.8", "outline_width": "0"}))
			rule_fill.setFilterExpression("\"ptype\" <> 'axis_center' AND \"cut_fill\" = 'fill'")
			rule_fill.setLabel("Auftrag (fill)")
			root.appendChild(rule_fill)
			# Achse rule for axis points — always show as white 1.1 pt dots
			rule_achse = QgsRuleBasedRenderer.Rule(QgsMarkerSymbol.createSimple({"color": "white", "outline_color": "white", "size": "1.1", "outline_width": "0", "name": "circle"}))
			rule_achse.setFilterExpression("\"ptype\" = 'axis_center'")
			rule_achse.setLabel("Achse (axis)")
			root.appendChild(rule_achse)
			layer.setRenderer(QgsRuleBasedRenderer(root))
			layer.triggerRepaint()

		def _style_line_layer(layer_id):
			layer = context.temporaryLayerStore().mapLayer(layer_id) or QgsProject.instance().mapLayer(layer_id)
			if not layer:
				return
			sym = QgsLineSymbol.createSimple({"color": "white", "width": "0.2"})
			layer.setRenderer(QgsSingleSymbolRenderer(sym))
			layer.triggerRepaint()


		_style_point_layer(points_id)
		try:
			feedback.pushInfo(f"Points added: {points_added}")
			layer_pts = context.temporaryLayerStore().mapLayer(points_id) or QgsProject.instance().mapLayer(points_id)
			if layer_pts:
				feedback.pushInfo(f"Points layer featureCount: {layer_pts.featureCount()}")
		except Exception:
			pass
		# Style lines
		for _lid in [segments_id, chain_c_id, transect_id]:
			_style_line_layer(_lid)

		# (Removed auto-queueing of layers; rely on standard Processing behavior)

		return {
			self.OUT_POINTS: points_id,
			self.OUT_BC_SEGMENTS: segments_id,
			self.OUT_AREAS: points_id,
			self.OUT_CHAIN_C: chain_c_id,
			self.OUT_CHAIN_EDGES: chain_c_id,
		}

	def name(self):
		# Unique (machine) name of this algorithm inside QGIS Processing.
		return "profiles_and_slope_intersections"

	def displayName(self):
		# Human-readable title shown in the toolbox.
		return "TRAILISM: Cross-Section Profiles"

	def shortHelpString(self):
		# One-paragraph help text shown in the Processing panel.
		return (
			"At each station along the input axis, the tool samples a DEM cross-section from s=-R (left) to s=+R (right) around the axis point (use to narrow down to a section of interest, ATTENTION: if R is too large, the tool may run out of memory, if R is too small, no intersection of sideslope with DEM is possible and C points in outputs will be missing!). The axis elevation z0 is taken from DEM plus optional dz, unless overridden by full cut/fill modes. The trail surface is flat between edges (no crossfall). From each trail edge, a side slope with angle beta_1 to horizontal is extended outward (up or down is chosen to hit the DEM first). Intersections define C points and the horizontal distance slope_run.\n"
			"Parameters: STEP_M (station spacing, m); WIDTH_FIELD or WIDTH_FALLBACK (trail width, m); PROFILE_RANGE_M R (m) for sampling extent; BETA1_DEG (deg) controls side slope tan(beta_1)=rise/run; DEM_OFFSET (m) shifts DEM vertically; AXIS_DZ (m) adjusts z0 as DEM+dz (ignored if full cut/fill); AXIS_FULL_CUT/FILL set z0 to min/max DEM at edges; SURF_* parameters only affect SVG visualization (trapezoid body h_A, side angle to horizontal, optional top crossfall q). SVG options allow plotting only selected stations with 1:1 scale. 'SVG: use all stations' will export all stations but is capped to a maximum of 100 SVGs per run.\n"
			"Outputs: points (axis center, edges, C points) with attributes including cut/fill kind, slope_run, z, z_base, z_diff, per-interval signed areas (area_cut>0, area_fill<0). Lines: BC segments per side; chains for C (left/right) and trail edges."
		)

	def createInstance(self):
		# Required for QGIS to instantiate the algorithm when loading the provider.
		return GenerateProfilesAndSlopeIntersections()