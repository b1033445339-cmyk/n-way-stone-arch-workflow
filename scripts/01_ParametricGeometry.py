"""Stage 01 - deterministic parametric geometry for N radial half arches.

Grasshopper Python inputs:
    N                  int          Number of half arches / supports (>= 3)
    H                  float        Overall crown top height above Z=0
    RotationAngles     list[float]  One plan angle in degrees per half arch
    Spans              list[float]  One crown-to-support span per half arch
    FootWidths         list[float]  One transverse width per half arch
    CrownDepth         float        Minimum vertical depth at the crown
    CrownFilletRadius  float        Plan fillet radius at crown intersections
    TopSlopePercent    float        Optional outward top fall, 0 to 5 percent
    CrownBlockRadius   float        Shared crown-zone / first-joint control

Use Item Access for N, H, CrownDepth, CrownFilletRadius, TopSlopePercent and
CrownBlockRadius. Use List Access for RotationAngles, Spans and FootWidths.
Empty RotationAngles gives an even 360-degree distribution. One Span or
FootWidth is broadcast to all arches. Unconnected crown inputs use their
DEFAULT values. An unconnected TopSlopePercent preserves the historical
horizontal top (0 percent); an unconnected CrownBlockRadius uses the same
1.60 m default as Stage 06.

Outputs:
    WholeSolid, CenterLines, SupportPoints, Valid, Report

This component performs geometry generation only. It does not use Kangaroo and
does not claim structural equilibrium. Every half arch is a deterministic
radial wall with a quadratic lower arch edge, an exact support location and a
constant transverse FootWidth. At zero slope all branches retain the historical
common horizontal top plane at Z=H. The wall is deepest at the support and has
the internal CrownDepth at the common crown. CrownFilletRadius rounds the
actual plan intersections of the half arches; it never adds a sphere or
ellipsoid. Block segmentation belongs to a later stage. A positive
TopSlopePercent retains a horizontal central cap of CrownBlockRadius and then
lowers each branch top radially outward at the requested gradient; it does not
cut the later crown block. The same CrownBlockRadius value is intended to feed
Stage 06 in parallel; no Stage 06 output is fed back into Stage 01.
"""

import math
import Rhino
import Rhino.Geometry as rg


DEFAULT_SPAN = 10.0
DEFAULT_FOOT_WIDTH = 1.0
DEFAULT_CROWN_DEPTH = 0.75
DEFAULT_CROWN_FILLET_RADIUS = 0.30
DEFAULT_TOP_SLOPE_PERCENT = 0.0
MAX_TOP_SLOPE_PERCENT = 5.0
DEFAULT_CROWN_BLOCK_RADIUS = 1.60


def _as_list(value):
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return [value]


def _document_tolerances():
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return 1.0e-3, math.radians(1.0)
    return (
        document.ModelAbsoluteTolerance,
        math.radians(document.ModelAngleToleranceDegrees),
    )


def _broadcast(values, count, default, name):
    prepared = _as_list(values)
    if not prepared:
        prepared = [default]
    if len(prepared) == 1:
        prepared *= count
    if len(prepared) != count:
        raise ValueError(
            "{0} must contain one value or exactly N values.".format(name)
        )
    prepared = [float(value) for value in prepared]
    if any(value <= 0.0 for value in prepared):
        raise ValueError("Every {0} value must be greater than zero.".format(
            name
        ))
    return prepared


def _top_slope_fraction(top_slope_percent):
    value = (
        DEFAULT_TOP_SLOPE_PERCENT
        if top_slope_percent is None
        else float(top_slope_percent)
    )
    if not math.isfinite(value):
        raise ValueError("TopSlopePercent must be finite.")
    if value < 0.0 or value > MAX_TOP_SLOPE_PERCENT:
        raise ValueError(
            "TopSlopePercent must be between 0 and {0:.6g}."
            .format(MAX_TOP_SLOPE_PERCENT)
        )
    return value / 100.0


def _crown_block_radius(crown_block_radius):
    value = (
        DEFAULT_CROWN_BLOCK_RADIUS
        if crown_block_radius is None
        else float(crown_block_radius)
    )
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("CrownBlockRadius must be a finite positive number.")
    return value


def _prepare_inputs(arch_count, height, rotation_angles, spans, foot_widths,
                    crown_depth, crown_fillet_radius):
    arch_count = int(arch_count)
    if arch_count < 3:
        raise ValueError("N must be at least 3.")

    height = float(height)
    crown_depth = (
        DEFAULT_CROWN_DEPTH
        if crown_depth is None
        else float(crown_depth)
    )
    if crown_depth <= 0.0:
        raise ValueError("CrownDepth must be greater than zero.")
    if height <= crown_depth:
        raise ValueError(
            "H must be greater than CrownDepth={0:.6g}."
            .format(crown_depth)
        )
    crown_fillet_radius = (
        DEFAULT_CROWN_FILLET_RADIUS
        if crown_fillet_radius is None
        else float(crown_fillet_radius)
    )
    if crown_fillet_radius < 0.0:
        raise ValueError("CrownFilletRadius cannot be negative.")

    angles = _as_list(rotation_angles)
    if not angles:
        angles = [
            360.0 * float(index) / float(arch_count)
            for index in range(arch_count)
        ]
    if len(angles) != arch_count:
        raise ValueError("RotationAngles must contain exactly N values.")
    angles = [float(angle) % 360.0 for angle in angles]

    rounded = [round(angle, 9) for angle in angles]
    if len(set(rounded)) != arch_count:
        raise ValueError("RotationAngles must be unique modulo 360 degrees.")

    span_values = _broadcast(
        spans, arch_count, DEFAULT_SPAN, "Spans"
    )
    width_values = _broadcast(
        foot_widths, arch_count, DEFAULT_FOOT_WIDTH, "FootWidths"
    )
    return (
        arch_count,
        height,
        angles,
        span_values,
        width_values,
        crown_depth,
        crown_fillet_radius,
    )


def _quadratic_curve(start, control, end):
    curve = rg.NurbsCurve.Create(False, 2, [start, control, end])
    if curve is None or not curve.IsValid:
        raise RuntimeError("Rhino could not create a quadratic arch curve.")
    return curve


def _arch_curves(angle_degrees, span, height, depth):
    angle = math.radians(angle_degrees)
    radial = rg.Vector3d(math.cos(angle), math.sin(angle), 0.0)
    transverse = rg.Vector3d(-radial.Y, radial.X, 0.0)

    bottom_start = rg.Point3d(
        radial.X * 0.0,
        radial.Y * 0.0,
        height - depth,
    )
    bottom_control = rg.Point3d(
        radial.X * (0.5 * span),
        radial.Y * (0.5 * span),
        height - depth,
    )
    bottom_end = rg.Point3d(
        radial.X * span,
        radial.Y * span,
        0.0,
    )
    bottom = _quadratic_curve(
        bottom_start, bottom_control, bottom_end
    )

    top_start = rg.Point3d(0.0, 0.0, height)
    top_end = rg.Point3d(
        radial.X * span,
        radial.Y * span,
        height,
    )
    top = rg.LineCurve(top_start, top_end)

    center_start = rg.Point3d(0.0, 0.0, height - 0.5 * depth)
    center_control = rg.Point3d(
        radial.X * (0.5 * span),
        radial.Y * (0.5 * span),
        height - 0.5 * depth,
    )
    center_end = rg.Point3d(
        radial.X * span,
        radial.Y * span,
        0.5 * height,
    )
    center = _quadratic_curve(
        center_start, center_control, center_end
    )
    return bottom, top, center, radial, transverse


def _top_elevation(height, slope, crown_block_radius, distance):
    return height - slope * max(
        0.0, distance - crown_block_radius
    )


def _sloped_top_curve(radial, span, height, slope, crown_block_radius):
    start = rg.Point3d(0.0, 0.0, height)
    transition = rg.Point3d(
        radial.X * crown_block_radius,
        radial.Y * crown_block_radius,
        height,
    )
    end = rg.Point3d(
        radial.X * span,
        radial.Y * span,
        _top_elevation(height, slope, crown_block_radius, span),
    )
    curve = rg.PolyCurve()
    if not curve.Append(rg.LineCurve(start, transition)):
        raise RuntimeError("Rhino could not append the horizontal crown cap.")
    if not curve.Append(rg.LineCurve(transition, end)):
        raise RuntimeError("Rhino could not append the sloped branch top.")
    if not curve.IsValid:
        raise RuntimeError("The sloped branch top is invalid.")
    return curve


def _sloped_center_line(radial, span, height, depth, slope,
                        crown_block_radius):
    points = []
    for index in range(65):
        distance = span * index / 64.0
        ratio = distance / span
        bottom_z = (height - depth) * (1.0 - ratio * ratio)
        top_z = _top_elevation(
            height, slope, crown_block_radius, distance
        )
        points.append(
            rg.Point3d(
                radial.X * distance,
                radial.Y * distance,
                0.5 * (top_z + bottom_z),
            )
        )
    curve = rg.PolylineCurve(points)
    if not curve.IsValid:
        raise RuntimeError("The sloped geometric centre line is invalid.")
    return curve


def _closed_profile(bottom, top):
    outer = rg.LineCurve(bottom.PointAtEnd, top.PointAtEnd)
    reversed_top = top.DuplicateCurve()
    if not reversed_top.Reverse():
        raise RuntimeError("Rhino could not reverse an arch top curve.")
    crown = rg.LineCurve(reversed_top.PointAtEnd, bottom.PointAtStart)

    profile = rg.PolyCurve()
    for segment in (bottom, outer, reversed_top, crown):
        if not profile.Append(segment):
            raise RuntimeError("Rhino could not join an arch profile.")
    if not profile.IsValid or not profile.IsClosed:
        raise RuntimeError("The generated arch profile is not closed.")
    return profile


def _extrude_profile(profile, transverse, width, tolerance):
    shifted = profile.DuplicateCurve()
    shifted.Transform(
        rg.Transform.Translation(transverse * (-0.5 * width))
    )
    surface = rg.Surface.CreateExtrusion(shifted, transverse * width)
    if surface is None:
        raise RuntimeError("Rhino could not extrude an arch profile.")
    solid = surface.ToBrep().CapPlanarHoles(tolerance)
    if solid is None or not solid.IsValid or not solid.IsSolid:
        raise RuntimeError("Rhino could not create a closed arch solid.")
    if solid.SolidOrientation == rg.BrepSolidOrientation.Inward:
        solid.Flip()
    if solid.SolidOrientation != rg.BrepSolidOrientation.Outward:
        raise RuntimeError("An arch solid has an invalid orientation.")
    return solid


def _valid_union(result):
    if result is None:
        return []
    return [
        brep for brep in result
        if (
            brep is not None
            and brep.IsValid
            and brep.IsSolid
        )
    ]


def _plan_rectangle(angle_degrees, width, start_distance, end_distance):
    angle = math.radians(angle_degrees)
    radial = rg.Vector3d(math.cos(angle), math.sin(angle), 0.0)
    transverse = rg.Vector3d(-radial.Y, radial.X, 0.0)
    half_width = 0.5 * width
    points = [
        rg.Point3d(
            start_distance * radial.X - half_width * transverse.X,
            start_distance * radial.Y - half_width * transverse.Y,
            0.0,
        ),
        rg.Point3d(
            end_distance * radial.X - half_width * transverse.X,
            end_distance * radial.Y - half_width * transverse.Y,
            0.0,
        ),
        rg.Point3d(
            end_distance * radial.X + half_width * transverse.X,
            end_distance * radial.Y + half_width * transverse.Y,
            0.0,
        ),
        rg.Point3d(
            start_distance * radial.X + half_width * transverse.X,
            start_distance * radial.Y + half_width * transverse.Y,
            0.0,
        ),
    ]
    points.append(points[0])
    return rg.PolylineCurve(points)


def _single_plan_union(curves, name):
    result = rg.Curve.CreateBooleanUnion(curves)
    if result is None or len(result) != 1:
        count = -1 if result is None else len(result)
        raise RuntimeError(
            "{0} produced {1} separate plan regions.".format(name, count)
        )
    curve = result[0]
    if not curve.IsValid or not curve.IsClosed:
        raise RuntimeError("{0} is not a valid closed curve.".format(name))
    return curve


def _rounded_full_plan(angles, widths, spans, fillet_radius, tolerance,
                       angle_tolerance):
    strips = [
        _plan_rectangle(angle, width, 0.0, span)
        for angle, width, span in zip(angles, widths, spans)
    ]
    sharp = _single_plan_union(strips, "The full arch footprint union")
    if fillet_radius <= tolerance:
        return sharp

    cleanup_length = 0.15 * min(widths)
    cleaned = sharp.DuplicateCurve()
    cleaned.RemoveShortSegments(cleanup_length)

    fully_filleted = rg.Curve.CreateFilletCornersCurve(
        cleaned,
        fillet_radius,
        tolerance,
        angle_tolerance,
    )
    if (
        fully_filleted is None
        or not fully_filleted.IsValid
        or not fully_filleted.IsClosed
    ):
        raise ValueError(
            "CrownFilletRadius={0:.6g} is too large for the current "
            "angles and FootWidths.".format(fillet_radius)
        )
    rounded = _single_plan_union(
        [cleaned, fully_filleted],
        "The selective concave fillet union",
    )
    curved_segments = sum(
        1 for segment in rounded.DuplicateSegments()
        if not segment.IsLinear(tolerance)
    )
    if curved_segments < len(angles):
        raise ValueError(
            "CrownFilletRadius created only {0} curved transitions for "
            "N={1}. Reduce the radius or revise the angle spacing."
            .format(curved_segments, len(angles))
        )
    return rounded


def _extrude_closed_plan_curve(plan_curve, base_z, extrusion_height):
    plan = plan_curve.DuplicateCurve()
    plan.Transform(rg.Transform.Translation(0.0, 0.0, base_z))
    extrusion = rg.Extrusion.Create(plan, extrusion_height, True)
    if extrusion is None:
        raise RuntimeError("Rhino could not extrude the crown plan.")
    solid = extrusion.ToBrep()
    if solid is None or not solid.IsValid or not solid.IsSolid:
        raise RuntimeError("The crown plan is not a valid closed solid.")
    if solid.SolidOrientation == rg.BrepSolidOrientation.Inward:
        solid.Flip()
    return solid


def _create_crown_geometry(angles, widths, spans, height, depth,
                           fillet_radius, tolerance, angle_tolerance):
    rounded_plan = _rounded_full_plan(
        angles,
        widths,
        spans,
        fillet_radius,
        tolerance,
        angle_tolerance,
    )
    overlap = max(20.0 * tolerance, 1.0e-5)
    connector = _extrude_closed_plan_curve(
        rounded_plan,
        height - depth - overlap,
        depth + 2.0 * overlap,
    )
    return connector


def _create_sloped_crown_geometry(angles, widths, spans, height, depth,
                                  fillet_radius, tolerance,
                                  angle_tolerance, crown_block_radius):
    rounded_plan = _rounded_full_plan(
        angles,
        widths,
        spans,
        fillet_radius,
        tolerance,
        angle_tolerance,
    )
    overlap = max(20.0 * tolerance, 1.0e-5)
    full_connector = _extrude_closed_plan_curve(
        rounded_plan,
        height - depth - overlap,
        depth + overlap,
    )
    cylinder_plane = rg.Plane(rg.Plane.WorldXY)
    cylinder_plane.Origin = rg.Point3d(
        0.0, 0.0, height - depth - overlap
    )
    cylinder = rg.Cylinder(
        rg.Circle(cylinder_plane, crown_block_radius),
        depth + overlap,
    ).ToBrep(True, True)
    results = _valid_union(
        rg.Brep.CreateBooleanIntersection(
            full_connector, cylinder, tolerance
        )
    )
    if len(results) != 1:
        raise RuntimeError(
            "The horizontal crown cap produced {0} solids."
            .format(len(results))
        )
    return results[0]


def _sequential_union(solids, tolerance):
    current = solids[0].DuplicateBrep()
    for solid in solids[1:]:
        result = _valid_union(
            rg.Brep.CreateBooleanUnion(
                [current, solid], tolerance
            )
        )
        if len(result) != 1:
            return None
        current = result[0]
    return current


def _trim_to_height(solid, height, extent, tolerance):
    clipping_box = rg.Box(
        rg.Plane.WorldXY,
        rg.Interval(-extent, extent),
        rg.Interval(-extent, extent),
        rg.Interval(-extent, height),
    ).ToBrep()
    result = _valid_union(
        rg.Brep.CreateBooleanIntersection(
            solid, clipping_box, tolerance
        )
    )
    return result[0] if len(result) == 1 else None


def create_parametric_geometry(arch_count, height, rotation_angles,
                               spans, foot_widths, crown_depth=None,
                               crown_fillet_radius=None,
                               top_slope_percent=None,
                               crown_block_radius=None):
    top_slope = _top_slope_fraction(top_slope_percent)
    crown_block_radius = _crown_block_radius(crown_block_radius)
    (
        arch_count,
        height,
        angles,
        span_values,
        width_values,
        crown_depth,
        crown_fillet_radius,
    ) = _prepare_inputs(
        arch_count,
        height,
        rotation_angles,
        spans,
        foot_widths,
        crown_depth,
        crown_fillet_radius,
    )
    if (
        top_slope > 0.0
        and min(span_values) <= crown_block_radius
    ):
        raise ValueError(
            "Every Span must exceed CrownBlockRadius."
        )
    tolerance, angle_tolerance = _document_tolerances()

    center_lines = []
    support_points = []
    solids = []
    drops = []
    minimum_sampled_thickness = float("inf")

    for angle, span, width in zip(
        angles, span_values, width_values
    ):
        bottom, horizontal_top, horizontal_center, radial, transverse = _arch_curves(
            angle, span, height, crown_depth
        )
        if top_slope > 0.0:
            top = _sloped_top_curve(
                radial, span, height, top_slope, crown_block_radius
            )
            center = _sloped_center_line(
                radial, span, height, crown_depth, top_slope,
                crown_block_radius
            )
            drop = top_slope * (
                span - crown_block_radius
            )
            drops.append(drop)
            for sample_index in range(65):
                distance = span * sample_index / 64.0
                ratio = distance / span
                bottom_z = (
                    (height - crown_depth) * (1.0 - ratio * ratio)
                )
                thickness = (
                    _top_elevation(
                        height, top_slope, crown_block_radius, distance
                    )
                    - bottom_z
                )
                minimum_sampled_thickness = min(
                    minimum_sampled_thickness, thickness
                )
        else:
            top = horizontal_top
            center = horizontal_center
        profile = _closed_profile(bottom, top)
        solid = _extrude_profile(
            profile, transverse, width, tolerance
        )
        center_lines.append(center)
        support_points.append(
            rg.Point3d(radial.X * span, radial.Y * span, 0.0)
        )
        solids.append(solid)

    if top_slope > 0.0:
        crown_connector = _create_sloped_crown_geometry(
            angles,
            width_values,
            span_values,
            height,
            crown_depth,
            crown_fillet_radius,
            tolerance,
            angle_tolerance,
            crown_block_radius,
        )
        union_results = _valid_union(
            rg.Brep.CreateBooleanUnion(
                [crown_connector] + solids, tolerance
            )
        )
        whole_solid = union_results[0] if len(union_results) == 1 else None
    else:
        crown_connector = _create_crown_geometry(
            angles,
            width_values,
            span_values,
            height,
            crown_depth,
            crown_fillet_radius,
            tolerance,
            angle_tolerance,
        )
        raw_union = _sequential_union(
            [crown_connector] + solids, tolerance
        )
        extent = max(span_values) + 3.0 * max(width_values)
        whole_solid = (
            _trim_to_height(raw_union, height, extent, tolerance)
            if raw_union is not None
            else None
        )

    height_error = 0.0
    for solid in solids:
        height_error = max(
            height_error,
            abs(solid.GetBoundingBox(True).Max.Z - height),
        )

    all_solids_valid = (
        len(solids) == arch_count
        and all(solid.IsValid and solid.IsSolid for solid in solids)
    )
    valid = (
        all_solids_valid
        and whole_solid is not None
        and height_error <= max(10.0 * tolerance, 1.0e-8)
        and (
            top_slope <= 0.0
            or minimum_sampled_thickness > tolerance
        )
    )
    if top_slope > 0.0:
        report = (
            f"N={arch_count}; geometry_mode=deterministic; H={height:.6g}; "
            "top_mode=horizontal_cap_then_radial_outward_fall; "
            f"top_slope_percent={100.0 * top_slope:.6g}; "
            f"crown_block_radius={crown_block_radius:.6g}; "
            f"top_drops_m={[round(value, 6) for value in drops]}; "
            f"minimum_sampled_thickness={minimum_sampled_thickness:.6g}; "
            f"crown_depth={crown_depth:.6g}; "
            f"crown_fillet_radius={crown_fillet_radius:.6g}; "
            "fillet_mode=selective_concave_union; "
            f"center_lines={len(center_lines)}; "
            f"supports={len(support_points)}; "
            f"height_error={height_error:.6g}; valid={valid}"
        )
    else:
        report = (
            f"N={arch_count}; geometry_mode=deterministic; H={height:.6g}; "
            "top_mode=common_horizontal; "
            f"crown_depth={crown_depth:.6g}; "
            f"crown_fillet_radius={crown_fillet_radius:.6g}; "
            "fillet_mode=selective_concave_union; "
            f"center_lines={len(center_lines)}; "
            f"supports={len(support_points)}; "
            f"height_error={height_error:.6g}; valid={valid}"
        )
    return (
        whole_solid,
        center_lines,
        support_points,
        valid,
        report,
    )


try:
    _top_slope_component_input = TopSlopePercent
except NameError:
    _top_slope_component_input = None

try:
    _crown_block_radius_component_input = CrownBlockRadius
except NameError:
    _crown_block_radius_component_input = None

try:
    _component_inputs = (
        N,
        H,
        RotationAngles,
        Spans,
        FootWidths,
        CrownDepth,
        CrownFilletRadius,
        _top_slope_component_input,
        _crown_block_radius_component_input,
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        WholeSolid,
        CenterLines,
        SupportPoints,
        Valid,
        Report,
    ) = create_parametric_geometry(*_component_inputs)
