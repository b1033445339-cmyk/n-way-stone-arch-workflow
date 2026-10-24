#! python 3
"""Stage 07D - compression-only partial-contact pressure redistribution.

Grasshopper Python inputs:
    JointProfiles      list[Curve]  JointProfiles from Stage 06 (List)
    JointPlanes        list[Plane]  JointPlanes from Stage 06 (List)
    JointIDs           list[str]    JointIDs from Stage 07A (List)
    CaseAngles         list[float]  CaseAngles from Stage 07B (List)
    CaseJointIDs       list[str]    CaseJointIDs from Stage 07B (List)
    ContactPoints      list[Point]  Flattened Stage 07B points (List)
    ContactResultants  list[Vector] Flattened Stage 07B resultants (List)
    KernUtilisations   list[float]  Flattened Stage 07C values (List)

Outputs:
    ContactAreaRatios, MinPressures, MaxPressures, PressureCentroids,
    PartialContactCaseKeys, Valid, Report

Full-contact joints retain the linear Stage 07C pressure field. For joints
outside the no-tension kern, an affine compression field truncated at zero is
solved over a triangulation of the actual joint profile. Its total force and
pressure centroid reproduce the Stage 07B normal resultant. Pressures are MPa.

This stage does not compare pressure with a stone strength and does not check
friction or dowel capacity.
"""

import math

import Rhino
import Rhino.Geometry as rg
from scipy.optimize import least_squares


BOUNDARY_SAMPLES = 192
MOMENT_EPSILON = 1.0e-12
KERN_TOLERANCE = 1.0e-6
MAX_SOLVER_EVALUATIONS = 500


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]


def _document_tolerance():
    document = Rhino.RhinoDoc.ActiveDoc
    return document.ModelAbsoluteTolerance if document is not None else 1.0e-3


def _coerce(value, geometry_type):
    if isinstance(value, geometry_type):
        return value
    if hasattr(value, "Value") and isinstance(value.Value, geometry_type):
        return value.Value
    return None


def _plane_coordinates(point, plane):
    vector = point - plane.Origin
    return (
        rg.Vector3d.Multiply(vector, plane.XAxis),
        rg.Vector3d.Multiply(vector, plane.YAxis),
    )


def _world_point(plane, x, y):
    return plane.Origin + plane.XAxis * x + plane.YAxis * y


def _sample_boundary(curve):
    try:
        success, polyline = curve.TryGetPolyline()
    except Exception:
        success, polyline = False, None
    if success and polyline is not None and polyline.Count >= 4:
        points = [rg.Point3d(point) for point in polyline]
        if points[0].DistanceTo(points[-1]) <= 1.0e-9:
            points.pop()
        if len(points) >= 3:
            return points
    parameters = curve.DivideByCount(BOUNDARY_SAMPLES, False)
    if parameters is None or len(parameters) < 3:
        raise RuntimeError("Could not sample a JointProfile boundary.")
    return [curve.PointAt(parameter) for parameter in parameters]


def _polygon_moments(points):
    area2 = sx6 = sy6 = 0.0
    sxx12 = syy12 = sxy24 = 0.0
    for index, point in enumerate(points):
        following = points[(index + 1) % len(points)]
        x0, y0 = point
        x1, y1 = following
        cross = x0 * y1 - x1 * y0
        area2 += cross
        sx6 += (x0 + x1) * cross
        sy6 += (y0 + y1) * cross
        sxx12 += (x0 * x0 + x0 * x1 + x1 * x1) * cross
        syy12 += (y0 * y0 + y0 * y1 + y1 * y1) * cross
        sxy24 += (
            2.0 * x0 * y0 + x0 * y1 + x1 * y0 + 2.0 * x1 * y1
        ) * cross
    area = area2 / 2.0
    if abs(area) <= MOMENT_EPSILON:
        raise RuntimeError("A polygon has zero area.")
    sign = 1.0 if area > 0.0 else -1.0
    return (
        area * sign,
        sx6 / 6.0 * sign,
        sy6 / 6.0 * sign,
        sxx12 / 12.0 * sign,
        syy12 / 12.0 * sign,
        sxy24 / 24.0 * sign,
    )


def _clean_polygon(points, tolerance):
    cleaned = []
    for point in points:
        if not cleaned or math.hypot(
            point[0] - cleaned[-1][0], point[1] - cleaned[-1][1]
        ) > tolerance:
            cleaned.append(point)
    if len(cleaned) > 1 and math.hypot(
        cleaned[0][0] - cleaned[-1][0],
        cleaned[0][1] - cleaned[-1][1],
    ) <= tolerance:
        cleaned.pop()
    changed = True
    while changed and len(cleaned) >= 3:
        changed = False
        reduced = []
        for index, point in enumerate(cleaned):
            previous = cleaned[(index - 1) % len(cleaned)]
            following = cleaned[(index + 1) % len(cleaned)]
            first_length = math.hypot(
                point[0] - previous[0], point[1] - previous[1]
            )
            second_length = math.hypot(
                following[0] - point[0], following[1] - point[1]
            )
            cross = abs(_cross(previous, point, following))
            if (
                first_length <= tolerance
                or second_length <= tolerance
                or cross <= tolerance * (first_length + second_length)
            ):
                changed = True
            else:
                reduced.append(point)
        if len(reduced) < 3:
            break
        cleaned = reduced
    return cleaned


def _cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _point_in_triangle(point, a, b, c, tolerance):
    return (
        _cross(a, b, point) >= -tolerance
        and _cross(b, c, point) >= -tolerance
        and _cross(c, a, point) >= -tolerance
    )


def _triangulate(boundary, tolerance):
    points = list(boundary)
    signed_area = 0.5 * sum(
        points[index][0] * points[(index + 1) % len(points)][1]
        - points[(index + 1) % len(points)][0] * points[index][1]
        for index in range(len(points))
    )
    if signed_area < 0.0:
        points.reverse()
    indices = list(range(len(points)))
    triangles = []
    guard = 0
    while len(indices) > 3 and guard < len(points) * len(points):
        guard += 1
        clipped = False
        for local_index, current in enumerate(indices):
            previous = indices[(local_index - 1) % len(indices)]
            following = indices[(local_index + 1) % len(indices)]
            a, b, c = points[previous], points[current], points[following]
            if _cross(a, b, c) <= tolerance:
                continue
            if any(
                _point_in_triangle(points[index], a, b, c, tolerance)
                for index in indices
                if index not in (previous, current, following)
            ):
                continue
            triangles.append([a, b, c])
            indices.pop(local_index)
            clipped = True
            break
        if not clipped:
            raise RuntimeError("Could not triangulate a JointProfile.")
    if len(indices) == 3:
        triangles.append([points[index] for index in indices])
    return triangles


def _clip_positive(polygon, ux, uy, cutoff):
    clipped = []
    for index, start in enumerate(polygon):
        end = polygon[(index + 1) % len(polygon)]
        start_value = ux * start[0] + uy * start[1] - cutoff
        end_value = ux * end[0] + uy * end[1] - cutoff
        start_inside = start_value >= -1.0e-12
        end_inside = end_value >= -1.0e-12
        if start_inside:
            clipped.append(start)
        if start_inside != end_inside:
            denominator = end_value - start_value
            if abs(denominator) > 1.0e-15:
                parameter = -start_value / denominator
                clipped.append((
                    start[0] + parameter * (end[0] - start[0]),
                    start[1] + parameter * (end[1] - start[1]),
                ))
    return clipped


def _positive_field_moments(triangles, ux, uy, cutoff):
    area = sx = sy = sxx = syy = sxy = 0.0
    for triangle in triangles:
        polygon = _clip_positive(triangle, ux, uy, cutoff)
        if len(polygon) < 3:
            continue
        values = _polygon_moments(polygon)
        area += values[0]
        sx += values[1]
        sy += values[2]
        sxx += values[3]
        syy += values[4]
        sxy += values[5]
    weighted = ux * sx + uy * sy - cutoff * area
    weighted_x = ux * sxx + uy * sxy - cutoff * sx
    weighted_y = ux * sxy + uy * syy - cutoff * sy
    if weighted <= MOMENT_EPSILON:
        return None
    return area, weighted, weighted_x / weighted, weighted_y / weighted


def _inverse_moment_action(x, y, second_x, second_y, product):
    determinant = second_x * second_y - product * product
    return (
        (second_y * x - product * y) / determinant,
        (-product * x + second_x * y) / determinant,
    )


def _prepare_joint_section(profile, plane, tolerance):
    world_boundary = _sample_boundary(profile)
    local_boundary = _clean_polygon(
        [_plane_coordinates(value, plane) for value in world_boundary],
        max(tolerance * 0.01, 1.0e-9),
    )
    area, sx, sy, sxx, syy, sxy = _polygon_moments(local_boundary)
    centroid = (sx / area, sy / area)
    centered_boundary = [
        (value[0] - centroid[0], value[1] - centroid[1])
        for value in local_boundary
    ]
    full_moments = _polygon_moments(centered_boundary)
    triangles = _triangulate(
        centered_boundary, max(tolerance * 1.0e-4, 1.0e-10)
    )
    normal = rg.Vector3d(plane.Normal)
    if not normal.Unitize():
        raise ValueError("JointPlane has an invalid normal.")
    return {
        "area": area,
        "centroid": centroid,
        "boundary": centered_boundary,
        "triangles": triangles,
        "full_moments": full_moments,
        "extent": max(
            max(math.hypot(*point) for point in centered_boundary),
            1.0e-6,
        ),
        "normal": normal,
    }


def _full_contact(boundary, normal_force, target, moments=None):
    area, sx, sy, sxx, syy, sxy = (
        moments if moments is not None else _polygon_moments(boundary)
    )
    centroid = (sx / area, sy / area)
    centered_target = (target[0] - centroid[0], target[1] - centroid[1])
    second_x = sxx - area * centroid[0] * centroid[0]
    second_y = syy - area * centroid[1] * centroid[1]
    product = sxy - area * centroid[0] * centroid[1]
    coefficient = _inverse_moment_action(
        centered_target[0], centered_target[1], second_x, second_y, product
    )
    terms = [
        1.0 / area
        + coefficient[0] * (point[0] - centroid[0])
        + coefficient[1] * (point[1] - centroid[1])
        for point in boundary
    ]
    return 1.0, normal_force * min(terms) / 1000.0, normal_force * max(terms) / 1000.0, target


def _partial_contact(
    boundary,
    triangles,
    normal_force,
    target,
    tolerance,
    area=None,
    extent=None,
    warm_start=None,
):
    if area is None:
        area = _polygon_moments(boundary)[0]
    if extent is None:
        extent = max(max(math.hypot(*point) for point in boundary), 1.0e-6)
    target_angle = math.atan2(target[1], target[0])
    best = None

    def state(parameters):
        theta, fraction = float(parameters[0]), float(parameters[1])
        ux, uy = math.cos(theta), math.sin(theta)
        projections = [ux * point[0] + uy * point[1] for point in boundary]
        minimum, maximum = min(projections), max(projections)
        cutoff = minimum + fraction * (maximum - minimum)
        moments = _positive_field_moments(triangles, ux, uy, cutoff)
        return ux, uy, cutoff, moments

    def residual(parameters):
        values = state(parameters)
        moments = values[3]
        if moments is None:
            return [1.0e3, 1.0e3]
        return [
            (moments[2] - target[0]) / extent,
            (moments[3] - target[1]) / extent,
        ]

    def solve(starts, current_best=None):
        best_result = current_best
        for start in starts:
            result = least_squares(
                residual,
                start,
                bounds=([-4.0 * math.pi, 0.0], [4.0 * math.pi, 0.999999]),
                xtol=1.0e-12,
                ftol=1.0e-12,
                gtol=1.0e-12,
                max_nfev=MAX_SOLVER_EVALUATIONS,
            )
            values = state(result.x)
            moments = values[3]
            if moments is None:
                continue
            error = math.hypot(
                moments[2] - target[0], moments[3] - target[1]
            )
            if best_result is None or error < best_result[0]:
                best_result = (error, values, moments, result)
        return best_result

    used_warm_start = warm_start is not None
    used_full_fallback = not used_warm_start
    if used_warm_start:
        warm_fraction = min(max(float(warm_start[1]), 0.0), 0.999999)
        warm_starts = [
            (float(warm_start[0]), warm_fraction),
            (target_angle, warm_fraction),
        ]
        best = solve(warm_starts)
        warm_tolerance = max(1.0e-9, tolerance * 1.0e-5)
        if best is None or best[0] > warm_tolerance:
            used_full_fallback = True

    if used_full_fallback:
        starts = []
        for angle_offset in (0.0, -0.35, 0.35, -0.75, 0.75, math.pi):
            for fraction in (0.10, 0.35, 0.60, 0.82):
                starts.append((target_angle + angle_offset, fraction))
        best = solve(starts, best)

    if best is None or best[0] > max(0.1 * tolerance, 1.0e-6):
        raise RuntimeError(
            "Partial-contact pressure centroid did not converge."
        )
    error, values, moments, result = best
    ux, uy, cutoff = values[:3]
    contact_area, weighted, centroid_x, centroid_y = moments
    scale = normal_force / weighted
    maximum_shape = max(
        max(0.0, ux * point[0] + uy * point[1] - cutoff)
        for point in boundary
    )
    return (
        contact_area / area,
        0.0,
        scale * maximum_shape / 1000.0,
        (centroid_x, centroid_y),
        error,
        int(result.nfev),
        (float(result.x[0]), float(result.x[1])),
        used_warm_start,
        used_full_fallback,
    )


def create_partial_contact_pressure(
    joint_profiles,
    joint_planes,
    contact_points,
    contact_resultants,
    kern_utilisations,
    joint_ids,
    _prepared_sections=None,
    _warm_starts=None,
    _solver_stats=None,
):
    profiles = [_coerce(value, rg.Curve) for value in _as_list(joint_profiles)]
    planes = [_coerce(value, rg.Plane) for value in _as_list(joint_planes)]
    points = [_coerce(value, rg.Point3d) for value in _as_list(contact_points)]
    resultants = []
    for value in _as_list(contact_resultants):
        vector = _coerce(value, rg.Vector3d)
        if vector is None:
            try:
                vector = rg.Vector3d(value)
            except Exception:
                vector = None
        resultants.append(vector)
    utilisations = [float(value) for value in _as_list(kern_utilisations)]
    identifiers = [str(value) for value in _as_list(joint_ids)]
    count = len(profiles)
    if not count or not all(
        len(values) == count
        for values in (planes, points, resultants, utilisations, identifiers)
    ):
        raise ValueError("All Stage 07D inputs must correspond one-to-one.")
    if any(value is None for value in profiles + planes + points + resultants):
        raise ValueError("Stage 07D contains invalid geometric input data.")
    if len(set(identifiers)) != count:
        raise ValueError("JointIDs must be unique.")
    if _prepared_sections is None:
        prepared_sections = [None] * count
    else:
        prepared_sections = list(_prepared_sections)
        if len(prepared_sections) != count:
            raise ValueError("Prepared Stage 07D sections must match JointIDs.")
    if _warm_starts is None:
        warm_starts = [None] * count
    else:
        warm_starts = _warm_starts
        if len(warm_starts) != count:
            raise ValueError("Stage 07D warm starts must match JointIDs.")
    solver_stats = _solver_stats if _solver_stats is not None else {}
    for key in ("partial_solves", "warm_attempts", "full_multistart_solves"):
        solver_stats.setdefault(key, 0)

    tolerance = _document_tolerance()
    ratios = []
    minimum_pressures = []
    maximum_pressures = []
    pressure_centroids = []
    partial_ids = []
    centroid_errors = []
    solver_evaluations = []
    failures = []

    values_by_joint = zip(
        profiles,
        planes,
        points,
        resultants,
        utilisations,
        prepared_sections,
    )
    for index, values in enumerate(values_by_joint):
        try:
            profile, plane, point, resultant, utilisation, section = values
            if section is None:
                section = _prepare_joint_section(profile, plane, tolerance)
            area = section["area"]
            centroid = section["centroid"]
            centered_boundary = section["boundary"]
            target_local = _plane_coordinates(point, plane)
            target = (
                target_local[0] - centroid[0],
                target_local[1] - centroid[1],
            )
            normal = section["normal"]
            normal_force = rg.Vector3d.Multiply(resultant, normal)
            if normal_force <= 0.0:
                raise ValueError("ContactResultant has no positive compression.")
            if utilisation <= 1.0 + KERN_TOLERANCE:
                ratio, minimum, maximum, reconstructed = _full_contact(
                    centered_boundary,
                    normal_force,
                    target,
                    section["full_moments"],
                )
                error = 0.0
                evaluations = 0
            else:
                partial_result = _partial_contact(
                    centered_boundary,
                    section["triangles"],
                    normal_force,
                    target,
                    tolerance,
                    area,
                    section["extent"],
                    warm_starts[index],
                )
                (
                    ratio,
                    minimum,
                    maximum,
                    reconstructed,
                    error,
                    evaluations,
                    warm_starts[index],
                    used_warm_start,
                    used_full_fallback,
                ) = partial_result
                solver_stats["partial_solves"] += 1
                solver_stats["warm_attempts"] += int(used_warm_start)
                solver_stats["full_multistart_solves"] += int(
                    used_full_fallback
                )
                partial_ids.append(identifiers[index])
            ratios.append(ratio)
            minimum_pressures.append(minimum)
            maximum_pressures.append(maximum)
            pressure_centroids.append(_world_point(
                plane,
                centroid[0] + reconstructed[0],
                centroid[1] + reconstructed[1],
            ))
            centroid_errors.append(error)
            solver_evaluations.append(evaluations)
        except Exception as error:
            failures.append((index, identifiers[index], str(error)))

    valid = (
        not failures
        and len(ratios) == count
        and all(0.0 < value <= 1.0 + 1.0e-9 for value in ratios)
        and all(value >= -1.0e-12 for value in minimum_pressures)
        and all(value > 0.0 for value in maximum_pressures)
    )
    report = (
        f"stage=partial_contact_pressure; joints={count}; "
        f"partial_contact_joints={len(partial_ids)}; "
        f"partial_contact_joint_ids={partial_ids}; "
        f"contact_area_ratio_min_max="
        f"({min(ratios) if ratios else 0.0:.6g},"
        f"{max(ratios) if ratios else 0.0:.6g}); "
        f"redistributed_pressure_min_max_MPa="
        f"({min(minimum_pressures) if minimum_pressures else 0.0:.6g},"
        f"{max(maximum_pressures) if maximum_pressures else 0.0:.6g}); "
        f"max_pressure_centroid_error_m="
        f"{max(centroid_errors) if centroid_errors else 0.0:.6g}; "
        f"max_solver_evaluations="
        f"{max(solver_evaluations) if solver_evaluations else 0}; "
        f"failures={failures}; valid={valid}; "
        "pressure_model=compression_only_affine_partial_contact; "
        "stone_strength_crushing_friction_and_pin_capacity_not_checked"
    )
    return (
        ratios,
        minimum_pressures,
        maximum_pressures,
        pressure_centroids,
        partial_ids,
        valid,
        report,
    )


def _case_key(case_index, angle, identifier):
    return "case[{0}]/angle[{1:.6g}]/{2}".format(
        case_index, angle, identifier
    )


def create_directional_partial_contact_pressure(
    joint_profiles,
    joint_planes,
    joint_ids,
    case_angles,
    case_joint_ids,
    contact_points,
    contact_resultants,
    kern_utilisations,
):
    identifiers = [str(value) for value in _as_list(joint_ids)]
    try:
        angles = [float(value) for value in _as_list(case_angles)]
        utilisations = [float(value) for value in _as_list(kern_utilisations)]
    except (TypeError, ValueError):
        raise ValueError(
            "CaseAngles and KernUtilisations must contain numeric values."
        )
    row_ids = [str(value) for value in _as_list(case_joint_ids)]
    points = _as_list(contact_points)
    resultants = _as_list(contact_resultants)
    joint_count = len(identifiers)
    expected_rows = len(angles) * joint_count
    if (
        not joint_count
        or not angles
        or len(set(identifiers)) != joint_count
        or any(
            len(values) != expected_rows
            for values in (row_ids, points, resultants, utilisations)
        )
    ):
        raise ValueError(
            "Stage 07D directional inputs must contain one complete "
            "JointIDs-sized block for every CaseAngles item."
        )

    profiles = [_coerce(value, rg.Curve) for value in _as_list(joint_profiles)]
    planes = [_coerce(value, rg.Plane) for value in _as_list(joint_planes)]
    if (
        len(profiles) != joint_count
        or len(planes) != joint_count
        or any(value is None for value in profiles + planes)
    ):
        raise ValueError(
            "JointProfiles, JointPlanes and JointIDs must correspond one-to-one."
        )
    tolerance = _document_tolerance()
    prepared_sections = []
    for profile, plane in zip(profiles, planes):
        try:
            prepared_sections.append(
                _prepare_joint_section(profile, plane, tolerance)
            )
        except Exception:
            prepared_sections.append(None)
    warm_starts = [None] * joint_count
    solver_stats = {}

    ratios = []
    minimum_pressures = []
    maximum_pressures = []
    pressure_centroids = []
    partial_keys = []
    failures = []
    for case_index, angle in enumerate(angles):
        start = case_index * joint_count
        end = start + joint_count
        if row_ids[start:end] != identifiers:
            failures.append((case_index, angle, "case_joint_id_order_mismatch"))
            continue
        result = create_partial_contact_pressure(
            profiles,
            planes,
            points[start:end],
            resultants[start:end],
            utilisations[start:end],
            identifiers,
            prepared_sections,
            warm_starts,
            solver_stats,
        )
        if not result[5]:
            failures.append((case_index, angle, result[6]))
            continue
        ratios.extend(result[0])
        minimum_pressures.extend(result[1])
        maximum_pressures.extend(result[2])
        pressure_centroids.extend(result[3])
        partial_keys.extend(
            _case_key(case_index, angle, identifier)
            for identifier in result[4]
        )

    valid = (
        not failures
        and len(ratios)
        == len(minimum_pressures)
        == len(maximum_pressures)
        == len(pressure_centroids)
        == expected_rows
        and all(0.0 < value <= 1.0 + 1.0e-9 for value in ratios)
        and all(value >= -1.0e-12 for value in minimum_pressures)
    )
    report = (
        f"stage=directional_partial_contact_pressure; "
        f"angle_cases={len(angles)}; joints_per_case={joint_count}; "
        f"flattened_rows={len(ratios)}; "
        f"partial_contact_case_rows={len(partial_keys)}; "
        f"contact_area_ratio_min_max="
        f"({min(ratios) if ratios else 0.0:.6g},"
        f"{max(ratios) if ratios else 0.0:.6g}); "
        f"redistributed_pressure_min_max_MPa="
        f"({min(minimum_pressures) if minimum_pressures else 0.0:.6g},"
        f"{max(maximum_pressures) if maximum_pressures else 0.0:.6g}); "
        f"partial_solves={solver_stats.get('partial_solves', 0)}; "
        f"warm_start_attempts={solver_stats.get('warm_attempts', 0)}; "
        f"full_multistart_solves="
        f"{solver_stats.get('full_multistart_solves', 0)}; "
        f"case_failures={failures}; valid={valid}; flattening=case_major; "
        "pressure_model=compression_only_affine_partial_contact"
    )
    return (
        ratios,
        minimum_pressures,
        maximum_pressures,
        pressure_centroids,
        partial_keys,
        valid,
        report,
    )


try:
    _component_inputs = (
        JointProfiles,
        JointPlanes,
        JointIDs,
        CaseAngles,
        CaseJointIDs,
        ContactPoints,
        ContactResultants,
        KernUtilisations,
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        ContactAreaRatios,
        MinPressures,
        MaxPressures,
        PressureCentroids,
        PartialContactCaseKeys,
        Valid,
        Report,
    ) = create_directional_partial_contact_pressure(*_component_inputs)
