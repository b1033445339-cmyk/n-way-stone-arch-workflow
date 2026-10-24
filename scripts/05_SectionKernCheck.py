#! python 3
"""Stage 05 - no-tension section-kern check for the TNA network.

Grasshopper Python inputs:
    WholeSolid  Brep        WholeSolid from 01_ParametricGeometry (Item)
    TNAEdges    list[Line]  TNAEdges from 04_TNAEnvelopeCheck (List)

Outputs:
    SectionKerns, CriticalPoints, Utilisations, ScreeningPass, Report, Advice

At the midpoint of every TNA edge, the actual solid is cut by a plane normal
to that edge. The planar section containing the thrust point is sampled and
its exact linear-elastic no-tension kern is constructed numerically. For a
rectangular section under one-axis eccentricity this reduces to the familiar
middle-third rule. In 3D it is intentionally not replaced by an arbitrary
pipe radius.

This is a pre-segmentation screening check, not a final structural pass/fail.
Utilisations correspond one-to-one with TNAEdges. Values above 1 identify
locations that should not be treated as unreinforced dry joints. Stone
strength, exact block weights, joint friction, local crushing, pins and
construction sequence are checked only after the cutting planes exist.
"""

import math

import Rhino
import Rhino.Geometry as rg


BOUNDARY_SAMPLES = 96
MOMENT_EPSILON = 1.0e-12
STRESS_TOLERANCE = 1.0e-7


def _as_list(value):
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return [value]


def _document_tolerance():
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return 1.0e-3
    return document.ModelAbsoluteTolerance


def _coerce_brep(value):
    if isinstance(value, rg.Brep):
        return value
    if hasattr(value, "Value") and isinstance(value.Value, rg.Brep):
        return value.Value
    return None


def _coerce_line(value):
    if isinstance(value, rg.Line):
        return rg.Line(value.From, value.To)
    if isinstance(value, rg.LineCurve):
        return rg.Line(value.PointAtStart, value.PointAtEnd)
    if hasattr(value, "Value"):
        return _coerce_line(value.Value)
    return None


def _point_key(point, tolerance):
    scale = 1.0 / max(tolerance, 1.0e-9)
    return (
        int(round(point.X * scale)),
        int(round(point.Y * scale)),
        int(round(point.Z * scale)),
    )


def _edge_branch_indices(lines, tolerance):
    key_to_index = {}
    indexed_edges = []
    for line in lines:
        indices = []
        for point in (line.From, line.To):
            key = _point_key(point, tolerance)
            if key not in key_to_index:
                key_to_index[key] = len(key_to_index)
            indices.append(key_to_index[key])
        indexed_edges.append(tuple(indices))

    adjacency = [[] for _ in range(len(key_to_index))]
    for edge_index, (a, b) in enumerate(indexed_edges):
        adjacency[a].append((b, edge_index))
        adjacency[b].append((a, edge_index))
    crown = max(
        range(len(adjacency)),
        key=lambda index: len(adjacency[index]),
    )
    edge_to_branch = {}
    for branch_index, (neighbour, first_edge) in enumerate(
        adjacency[crown]
    ):
        previous = crown
        current = neighbour
        edge_to_branch[first_edge] = branch_index
        while current != crown and len(adjacency[current]) == 2:
            next_items = [
                item for item in adjacency[current]
                if item[0] != previous
            ]
            if not next_items:
                break
            next_node, next_edge = next_items[0]
            edge_to_branch[next_edge] = branch_index
            previous, current = current, next_node
    return [
        edge_to_branch.get(edge_index, -1)
        for edge_index in range(len(lines))
    ]


def _section_plane(line):
    midpoint = line.PointAt(0.5)
    normal = rg.Vector3d(line.Direction)
    if not normal.Unitize():
        raise RuntimeError("A TNA edge has zero length.")
    reference = rg.Vector3d.ZAxis
    if abs(rg.Vector3d.Multiply(normal, reference)) > 0.92:
        reference = rg.Vector3d.XAxis
    x_axis = rg.Vector3d.CrossProduct(reference, normal)
    if not x_axis.Unitize():
        raise RuntimeError("Could not construct a section plane.")
    y_axis = rg.Vector3d.CrossProduct(normal, x_axis)
    if not y_axis.Unitize():
        raise RuntimeError("Could not construct a section plane.")
    return midpoint, rg.Plane(midpoint, x_axis, y_axis)


def _joined_section_curves(solid, plane, tolerance):
    result = rg.Intersect.Intersection.BrepPlane(
        solid, plane, tolerance
    )
    if not result or not result[0]:
        return []
    raw_curves = list(result[1]) if result[1] else []
    joined = list(rg.Curve.JoinCurves(raw_curves, tolerance))
    return [
        curve for curve in joined
        if curve is not None and curve.IsValid and curve.IsClosed
    ]


def _curve_area(curve):
    properties = rg.AreaMassProperties.Compute(curve)
    if properties is None:
        return 0.0
    area = abs(properties.Area)
    properties.Dispose()
    return area


def _containing_section(curves, point, plane, tolerance):
    containing = []
    for curve in curves:
        status = curve.Contains(point, plane, tolerance)
        if status in (
            rg.PointContainment.Inside,
            rg.PointContainment.Coincident,
        ):
            containing.append(curve)
    if not containing:
        for curve in curves:
            curve.Dispose()
        raise RuntimeError(
            "No section containing the thrust point was found."
        )
    section = max(containing, key=_curve_area)
    for curve in curves:
        if curve is not section:
            curve.Dispose()
    return section


def _plane_coordinates(point, plane):
    vector = point - plane.Origin
    return (
        rg.Vector3d.Multiply(vector, plane.XAxis),
        rg.Vector3d.Multiply(vector, plane.YAxis),
    )


def _world_point(plane, x, y):
    return (
        plane.Origin
        + plane.XAxis * x
        + plane.YAxis * y
    )


def _sample_closed_curve(curve):
    parameters = curve.DivideByCount(BOUNDARY_SAMPLES, False)
    if not parameters:
        return []
    return [curve.PointAt(parameter) for parameter in parameters]


def _signed_area_and_centroid(points):
    twice_area = 0.0
    centroid_x = 0.0
    centroid_y = 0.0
    for index, point in enumerate(points):
        next_point = points[(index + 1) % len(points)]
        cross = (
            point[0] * next_point[1]
            - next_point[0] * point[1]
        )
        twice_area += cross
        centroid_x += (point[0] + next_point[0]) * cross
        centroid_y += (point[1] + next_point[1]) * cross
    area = 0.5 * twice_area
    if abs(area) <= MOMENT_EPSILON:
        raise RuntimeError("A sampled section has zero area.")
    centroid_x /= 6.0 * area
    centroid_y /= 6.0 * area
    return area, (centroid_x, centroid_y)


def _section_moments(points):
    area, centroid = _signed_area_and_centroid(points)
    if area < 0.0:
        points = list(reversed(points))
        area, centroid = _signed_area_and_centroid(points)
    centered = [
        (point[0] - centroid[0], point[1] - centroid[1])
        for point in points
    ]
    second_x = 0.0
    second_y = 0.0
    product = 0.0
    for index, point in enumerate(centered):
        next_point = centered[(index + 1) % len(centered)]
        cross = (
            point[0] * next_point[1]
            - next_point[0] * point[1]
        )
        second_x += (
            point[0] * point[0]
            + point[0] * next_point[0]
            + next_point[0] * next_point[0]
        ) * cross
        second_y += (
            point[1] * point[1]
            + point[1] * next_point[1]
            + next_point[1] * next_point[1]
        ) * cross
        product += (
            2.0 * point[0] * point[1]
            + point[0] * next_point[1]
            + next_point[0] * point[1]
            + 2.0 * next_point[0] * next_point[1]
        ) * cross
    second_x /= 12.0
    second_y /= 12.0
    product /= 24.0
    determinant = second_x * second_y - product * product
    if determinant <= MOMENT_EPSILON:
        raise RuntimeError("A sampled section has singular area moments.")
    return area, centroid, centered, second_x, second_y, product


def _inverse_moment_action(x, y, second_x, second_y, product):
    determinant = second_x * second_y - product * product
    return (
        (second_y * x - product * y) / determinant,
        (-product * x + second_x * y) / determinant,
    )


def _clip_half_plane(polygon, normal_x, normal_y, limit):
    if not polygon:
        return []
    clipped = []
    for index, start in enumerate(polygon):
        end = polygon[(index + 1) % len(polygon)]
        start_value = normal_x * start[0] + normal_y * start[1]
        end_value = normal_x * end[0] + normal_y * end[1]
        start_inside = start_value >= limit - 1.0e-12
        end_inside = end_value >= limit - 1.0e-12
        if start_inside:
            clipped.append(start)
        if start_inside != end_inside:
            denominator = end_value - start_value
            if abs(denominator) > 1.0e-15:
                factor = (limit - start_value) / denominator
                clipped.append(
                    (
                        start[0] + factor * (end[0] - start[0]),
                        start[1] + factor * (end[1] - start[1]),
                    )
                )
    return clipped


def _clean_polygon(polygon, tolerance):
    """Remove repeated and collinear vertices created by repeated clipping."""
    if not polygon:
        return []
    cleaned = []
    for point in polygon:
        if not cleaned:
            cleaned.append(point)
            continue
        distance = math.hypot(
            point[0] - cleaned[-1][0],
            point[1] - cleaned[-1][1],
        )
        if distance > tolerance:
            cleaned.append(point)
    if (
        len(cleaned) > 1
        and math.hypot(
            cleaned[0][0] - cleaned[-1][0],
            cleaned[0][1] - cleaned[-1][1],
        ) <= tolerance
    ):
        cleaned.pop()

    changed = True
    while changed and len(cleaned) >= 3:
        changed = False
        reduced = []
        count = len(cleaned)
        for index, point in enumerate(cleaned):
            previous = cleaned[(index - 1) % count]
            following = cleaned[(index + 1) % count]
            first_x = point[0] - previous[0]
            first_y = point[1] - previous[1]
            second_x = following[0] - point[0]
            second_y = following[1] - point[1]
            first_length = math.hypot(first_x, first_y)
            second_length = math.hypot(second_x, second_y)
            if first_length <= tolerance or second_length <= tolerance:
                changed = True
                continue
            cross = abs(first_x * second_y - first_y * second_x)
            if cross <= tolerance * (first_length + second_length):
                changed = True
                continue
            reduced.append(point)
        if len(reduced) < 3:
            break
        cleaned = reduced
    return cleaned


def _kern_polygon(area, boundary, second_x, second_y, product):
    extent = max(
        max(math.hypot(point[0], point[1]) for point in boundary),
        1.0e-6,
    )
    polygon = [
        (-2.0 * extent, -2.0 * extent),
        (2.0 * extent, -2.0 * extent),
        (2.0 * extent, 2.0 * extent),
        (-2.0 * extent, 2.0 * extent),
    ]
    limit = -1.0 / area
    cleanup_tolerance = max(extent * 1.0e-9, 1.0e-12)
    for point in boundary:
        normal_x, normal_y = _inverse_moment_action(
            point[0],
            point[1],
            second_x,
            second_y,
            product,
        )
        polygon = _clip_half_plane(
            polygon, normal_x, normal_y, limit
        )
        polygon = _clean_polygon(polygon, cleanup_tolerance)
        if len(polygon) < 3:
            return []
    return _clean_polygon(polygon, cleanup_tolerance)


def _section_check(solid, line, tolerance):
    thrust_point, plane = _section_plane(line)
    section_curves = _joined_section_curves(
        solid, plane, tolerance
    )
    section = _containing_section(
        section_curves, thrust_point, plane, tolerance
    )
    world_boundary = _sample_closed_curve(section)
    if len(world_boundary) < 8:
        raise RuntimeError("The solid section could not be sampled.")
    local_boundary = [
        _plane_coordinates(point, plane)
        for point in world_boundary
    ]
    (
        area,
        centroid,
        centered_boundary,
        second_x,
        second_y,
        product,
    ) = _section_moments(local_boundary)
    thrust_local = _plane_coordinates(thrust_point, plane)
    eccentricity = (
        thrust_local[0] - centroid[0],
        thrust_local[1] - centroid[1],
    )
    coefficient = _inverse_moment_action(
        eccentricity[0],
        eccentricity[1],
        second_x,
        second_y,
        product,
    )
    minimum_term = min(
        coefficient[0] * point[0] + coefficient[1] * point[1]
        for point in centered_boundary
    )
    utilisation = max(0.0, -area * minimum_term)
    passed = utilisation <= 1.0 + STRESS_TOLERANCE
    kern = _kern_polygon(
        area,
        centered_boundary,
        second_x,
        second_y,
        product,
    )
    if len(kern) < 3:
        section.Dispose()
        raise RuntimeError(
            "The calculated section kern polygon is empty."
        )
    kern_points = [
        _world_point(
            plane,
            centroid[0] + point[0],
            centroid[1] + point[1],
        )
        for point in kern
    ]
    kern_points.append(kern_points[0])
    kern_curve = rg.PolylineCurve(kern_points)
    if not kern_curve.IsValid or not kern_curve.IsClosed:
        section.Dispose()
        raise RuntimeError(
            "The calculated section kern curve is invalid."
        )
    section.Dispose()
    return passed, utilisation, kern_curve, thrust_point, area


def create_section_kern_check(whole_solid, tna_edges):
    solid = _coerce_brep(whole_solid)
    lines = []
    for edge_index, value in enumerate(_as_list(tna_edges)):
        line = _coerce_line(value)
        if line is None or not line.IsValid or line.Length <= 0.0:
            raise ValueError(
                "TNAEdges contains an invalid line at index {0}."
                .format(edge_index)
            )
        lines.append(line)
    if solid is None or not solid.IsValid or not solid.IsSolid:
        raise ValueError("WholeSolid must be one valid closed Brep.")
    if not lines:
        raise ValueError("TNAEdges must contain valid lines.")

    tolerance = _document_tolerance()
    branches = _edge_branch_indices(
        lines, max(tolerance * 0.1, 1.0e-7)
    )
    safety_zones = []
    critical_points = []
    failures = []
    missing_sections = []
    utilisations = [None for line in lines]
    areas = []

    for edge_index, line in enumerate(lines):
        try:
            (
                passed,
                utilisation,
                kern_curve,
                thrust_point,
                area,
            ) = _section_check(solid, line, tolerance)
            if kern_curve is not None:
                safety_zones.append(kern_curve)
            utilisations[edge_index] = utilisation
            areas.append(area)
            if not passed:
                critical_points.append(thrust_point)
                failures.append(
                    (edge_index, branches[edge_index], utilisation)
                )
        except Exception as error:
            missing_sections.append((edge_index, str(error)))
            critical_points.append(line.PointAt(0.5))

    screening_pass = (
        len(safety_zones) == len(lines)
        and not failures
        and not missing_sections
    )
    evaluated_utilisations = [
        value for value in utilisations if value is not None
    ]
    worst_utilisation = (
        max(evaluated_utilisations)
        if evaluated_utilisations
        else float("inf")
    )
    branch_worst = {}
    for edge_index, branch_index, utilisation in failures:
        branch_worst[branch_index] = max(
            branch_worst.get(branch_index, 0.0), utilisation
        )
    ranked_branches = [
        branch + 1 for branch, value in sorted(
            branch_worst.items(),
            key=lambda item: item[1],
            reverse=True,
        )
        if branch >= 0
    ]

    if screening_pass:
        advice = (
            "PRE-SEGMENTATION PASS: every sampled baseline thrust "
            "resultant lies inside the actual section kern. Joints may be "
            "placed near these stations, but actual-joint pressure and sliding "
            "require separate assessment; block equilibrium is not established."
        )
    elif missing_sections:
        advice = (
            "PRE-SEGMENTATION DATA ERROR: section construction failed at "
            "{0} stations. Confirm that "
            "WholeSolid is one closed union and TNAEdges comes directly "
            "from Stage 04.".format(len(missing_sections))
        )
    else:
        advice = (
            "PRE-SEGMENTATION WARNING: {0} sampled stations lie outside "
            "the no-tension kern; affected arch indices={1}. This does "
            "not reject all networks or force states. The failed screen concerns "
            "this selected network under the linear full-section pressure assumption. "
            "Stage 06 selects cutting planes; partial-contact and sliding-deficit "
            "assessment may proceed without claiming to repair this screen failure.".format(
                len(failures), ranked_branches
            )
        )

    report = (
        "stage=pre_segmentation_screening; "
        "criterion=actual_section_no_tension_kern; "
        "rectangular_uniaxial_equivalent=middle_third; edges={0}; "
        "sections_evaluated={1}; safety_zones={2}; failed_sections={3}; "
        "missing_sections={4}; worst_utilisation={5:.6g}; "
        "limit=1; min_section_area={6:.6g}; max_section_area={7:.6g}; "
        "affected_arches={8}; screening_pass={9}; "
        "final_structural_validity=not_evaluated"
    ).format(
        len(lines),
        len(evaluated_utilisations),
        len(safety_zones),
        len(failures),
        len(missing_sections),
        worst_utilisation,
        min(areas) if areas else 0.0,
        max(areas) if areas else 0.0,
        ranked_branches,
        screening_pass,
    )
    return (
        safety_zones,
        critical_points,
        utilisations,
        screening_pass,
        report,
        advice,
    )


try:
    _component_inputs = (WholeSolid, TNAEdges)
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        SectionKerns,
        CriticalPoints,
        Utilisations,
        ScreeningPass,
        Report,
        Advice,
    ) = (
        create_section_kern_check(*_component_inputs)
    )
