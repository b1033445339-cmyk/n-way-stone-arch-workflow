#! python 3
"""Stage 07C - actual-joint eccentricity, opening and contact pressure.

Grasshopper Python inputs:
    JointProfiles      list[Curve]  JointProfiles from Stage 06 (List)
    JointPlanes        list[Plane]  JointPlanes from Stage 06 (List)
    JointIDs           list[str]    JointIDs from Stage 07A (List)
    CaseAngles         list[float]  CaseAngles from Stage 07B (List)
    CaseJointIDs       list[str]    CaseJointIDs from Stage 07B (List)
    ContactPoints      list[Point]  Flattened Stage 07B points (List)
    ContactResultants  list[Vector] Flattened Stage 07B resultants (List)

Outputs:
    JointKerns, KernUtilisations, MinPressures, MaxPressures,
    OpeningCaseKeys, NoOpeningPass, Valid, Report

The actual Stage 06 joint profile is treated as a planar section. A linear
normal-pressure field is reconstructed from the normal resultant and its
eccentricity. The exact linear-elastic no-tension kern is generated numerically
from the sampled profile. Pressures are returned in MPa.

Negative MinPressures indicate loss of full contact. In that case the reported
linear pressure range is a diagnostic full-contact field, not a nonlinear
partial-contact solution. Stone crushing, friction and dowel capacity are not
checked in this stage.
"""

import math

import Rhino
import Rhino.Geometry as rg


BOUNDARY_SAMPLES = 192
MOMENT_EPSILON = 1.0e-12
KERN_TOLERANCE = 1.0e-6


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
    if document is None:
        return 1.0e-3
    return document.ModelAbsoluteTolerance


def _coerce_curve(value):
    if isinstance(value, rg.Curve):
        return value
    if hasattr(value, "Value") and isinstance(value.Value, rg.Curve):
        return value.Value
    return None


def _coerce_plane(value):
    if isinstance(value, rg.Plane):
        return rg.Plane(value)
    if hasattr(value, "Value") and isinstance(value.Value, rg.Plane):
        return rg.Plane(value.Value)
    return None


def _coerce_point(value):
    if isinstance(value, rg.Point3d):
        return rg.Point3d(value)
    if hasattr(value, "Value") and isinstance(value.Value, rg.Point3d):
        return rg.Point3d(value.Value)
    return None


def _coerce_vector(value):
    if isinstance(value, rg.Vector3d):
        return rg.Vector3d(value)
    if hasattr(value, "Value") and isinstance(value.Value, rg.Vector3d):
        return rg.Vector3d(value.Value)
    try:
        return rg.Vector3d(value)
    except Exception:
        return None


def _plane_coordinates(point, plane):
    vector = point - plane.Origin
    return (
        rg.Vector3d.Multiply(vector, plane.XAxis),
        rg.Vector3d.Multiply(vector, plane.YAxis),
    )


def _world_point(plane, x, y):
    return plane.Origin + plane.XAxis * x + plane.YAxis * y


def _sample_closed_curve(curve):
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


def _signed_area_and_centroid(points):
    twice_area = 0.0
    centroid_x = 0.0
    centroid_y = 0.0
    for index, point in enumerate(points):
        following = points[(index + 1) % len(points)]
        cross = point[0] * following[1] - following[0] * point[1]
        twice_area += cross
        centroid_x += (point[0] + following[0]) * cross
        centroid_y += (point[1] + following[1]) * cross
    area = 0.5 * twice_area
    if abs(area) <= MOMENT_EPSILON:
        raise RuntimeError("A sampled JointProfile has zero area.")
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
        following = centered[(index + 1) % len(centered)]
        cross = point[0] * following[1] - following[0] * point[1]
        second_x += (
            point[0] * point[0]
            + point[0] * following[0]
            + following[0] * following[0]
        ) * cross
        second_y += (
            point[1] * point[1]
            + point[1] * following[1]
            + following[1] * following[1]
        ) * cross
        product += (
            2.0 * point[0] * point[1]
            + point[0] * following[1]
            + following[0] * point[1]
            + 2.0 * following[0] * following[1]
        ) * cross
    second_x /= 12.0
    second_y /= 12.0
    product /= 24.0
    determinant = second_x * second_y - product * product
    if determinant <= MOMENT_EPSILON:
        raise RuntimeError("A JointProfile has singular area moments.")
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
        start_inside = start_value <= limit + 1.0e-12
        end_inside = end_value <= limit + 1.0e-12
        if start_inside:
            clipped.append(start)
        if start_inside != end_inside:
            denominator = end_value - start_value
            if abs(denominator) > 1.0e-15:
                parameter = (limit - start_value) / denominator
                clipped.append(
                    (
                        start[0] + (end[0] - start[0]) * parameter,
                        start[1] + (end[1] - start[1]) * parameter,
                    )
                )
    return clipped


def _clean_polygon(polygon, tolerance):
    if not polygon:
        return []
    cleaned = []
    for point in polygon:
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
        result = []
        count = len(cleaned)
        for index, point in enumerate(cleaned):
            previous = cleaned[(index - 1) % count]
            following = cleaned[(index + 1) % count]
            cross = (
                (point[0] - previous[0]) * (following[1] - point[1])
                - (point[1] - previous[1]) * (following[0] - point[0])
            )
            scale = max(
                math.hypot(point[0] - previous[0], point[1] - previous[1]),
                math.hypot(following[0] - point[0], following[1] - point[1]),
                1.0,
            )
            if abs(cross) <= tolerance * scale:
                changed = True
            else:
                result.append(point)
        cleaned = result
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
    cleanup_tolerance = max(extent * 1.0e-9, 1.0e-12)
    for point in boundary:
        normal_x, normal_y = _inverse_moment_action(
            point[0], point[1], second_x, second_y, product
        )
        polygon = _clip_half_plane(
            polygon, -normal_x, -normal_y, 1.0 / area
        )
        polygon = _clean_polygon(polygon, cleanup_tolerance)
        if len(polygon) < 3:
            return []
    return _clean_polygon(polygon, cleanup_tolerance)


def _prepare_joint_section(profile, plane):
    world_boundary = _sample_closed_curve(profile)
    local_boundary = [
        _plane_coordinates(point, plane) for point in world_boundary
    ]
    (
        area,
        centroid,
        centered_boundary,
        second_x,
        second_y,
        product,
    ) = _section_moments(local_boundary)
    kern = _kern_polygon(
        area, centered_boundary, second_x, second_y, product
    )
    if len(kern) < 3:
        raise RuntimeError("The calculated actual-joint kern is empty.")
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
        raise RuntimeError("The calculated actual-joint kern is invalid.")
    normal = rg.Vector3d(plane.Normal)
    if not normal.Unitize():
        raise ValueError("A JointPlane has an invalid normal.")
    return {
        "area": area,
        "centroid": centroid,
        "centered_boundary": centered_boundary,
        "second_x": second_x,
        "second_y": second_y,
        "product": product,
        "kern": kern_curve,
        "normal": normal,
    }


def _joint_check(
    profile,
    plane,
    contact_point,
    resultant,
    tolerance,
    prepared_section=None,
):
    containment = profile.Contains(contact_point, plane, tolerance)
    if containment not in (
        rg.PointContainment.Inside,
        rg.PointContainment.Coincident,
    ):
        raise ValueError("A ContactPoint lies outside its JointProfile.")

    section = (
        prepared_section
        if prepared_section is not None
        else _prepare_joint_section(profile, plane)
    )
    area = section["area"]
    centroid = section["centroid"]
    centered_boundary = section["centered_boundary"]
    second_x = section["second_x"]
    second_y = section["second_y"]
    product = section["product"]
    contact_local = _plane_coordinates(contact_point, plane)
    eccentricity = (
        contact_local[0] - centroid[0],
        contact_local[1] - centroid[1],
    )
    coefficient = _inverse_moment_action(
        eccentricity[0],
        eccentricity[1],
        second_x,
        second_y,
        product,
    )
    pressure_terms = [
        1.0 / area + coefficient[0] * point[0] + coefficient[1] * point[1]
        for point in centered_boundary
    ]
    minimum_term = min(pressure_terms)
    maximum_term = max(pressure_terms)
    utilisation = max(
        0.0,
        -area * min(
            coefficient[0] * point[0] + coefficient[1] * point[1]
            for point in centered_boundary
        ),
    )
    no_opening = utilisation <= 1.0 + KERN_TOLERANCE

    normal = section["normal"]
    normal_force = rg.Vector3d.Multiply(resultant, normal)
    tangential = resultant - normal * normal_force
    if normal_force <= 0.0:
        raise ValueError("A ContactResultant has no positive compression.")
    minimum_pressure = normal_force * minimum_term / 1000.0
    maximum_pressure = normal_force * maximum_term / 1000.0

    return {
        "kern": section["kern"],
        "utilisation": utilisation,
        "no_opening": no_opening,
        "minimum_pressure": minimum_pressure,
        "maximum_pressure": maximum_pressure,
        "area": area,
        "eccentricity": math.hypot(*eccentricity),
        "normal_force": normal_force,
        "tangential_force": tangential.Length,
    }


def create_joint_opening_and_pressure(
    joint_profiles,
    joint_planes,
    contact_points,
    contact_resultants,
    joint_ids,
    _prepared_sections=None,
):
    profiles = []
    for index, value in enumerate(_as_list(joint_profiles)):
        profile = _coerce_curve(value)
        if profile is None or not profile.IsValid or not profile.IsClosed:
            raise ValueError(
                "JointProfiles contains invalid data at index {0}."
                .format(index)
            )
        profiles.append(profile)

    planes = []
    for index, value in enumerate(_as_list(joint_planes)):
        plane = _coerce_plane(value)
        if plane is None or not plane.IsValid:
            raise ValueError(
                "JointPlanes contains invalid data at index {0}."
                .format(index)
            )
        planes.append(plane)

    points = []
    for index, value in enumerate(_as_list(contact_points)):
        point = _coerce_point(value)
        if point is None or not point.IsValid:
            raise ValueError(
                "ContactPoints contains invalid data at index {0}."
                .format(index)
            )
        points.append(point)

    resultants = []
    for index, value in enumerate(_as_list(contact_resultants)):
        resultant = _coerce_vector(value)
        if resultant is None or not all(
            math.isfinite(component)
            for component in (resultant.X, resultant.Y, resultant.Z)
        ):
            raise ValueError(
                "ContactResultants contains invalid data at index {0}."
                .format(index)
            )
        resultants.append(resultant)

    identifiers = [str(value) for value in _as_list(joint_ids)]
    count = len(profiles)
    if not (
        len(planes) == count
        and len(points) == count
        and len(resultants) == count
        and len(identifiers) == count
    ):
        raise ValueError(
            "All Stage 07C inputs must correspond one-to-one with joints."
        )
    if not identifiers or len(set(identifiers)) != count:
        raise ValueError("JointIDs must contain unique identifiers.")
    if _prepared_sections is None:
        prepared_sections = [None] * count
    else:
        prepared_sections = list(_prepared_sections)
        if len(prepared_sections) != count:
            raise ValueError("Prepared Stage 07C sections must match JointIDs.")

    tolerance = _document_tolerance()
    records = []
    failures = []
    values_by_joint = zip(
        profiles,
        planes,
        points,
        resultants,
        prepared_sections,
    )
    for index, values in enumerate(values_by_joint):
        try:
            profile, plane, point, resultant, prepared = values
            records.append(_joint_check(
                profile,
                plane,
                point,
                resultant,
                tolerance,
                prepared,
            ))
        except Exception as error:
            failures.append((index, identifiers[index], str(error)))

    joint_kerns = [record["kern"] for record in records]
    utilisations = [record["utilisation"] for record in records]
    minimum_pressures = [record["minimum_pressure"] for record in records]
    maximum_pressures = [record["maximum_pressure"] for record in records]
    opening_ids = [
        identifiers[index]
        for index, record in enumerate(records)
        if not record["no_opening"]
    ] if not failures else []
    no_opening_pass = not failures and not opening_ids
    valid = not failures and len(records) == count

    report = (
        f"stage=joint_opening_and_pressure; joints={count}; "
        f"sections_evaluated={len(records)}; kerns={len(joint_kerns)}; "
        f"opening_joints={len(opening_ids)}; opening_joint_ids={opening_ids}; "
        f"no_opening_pass={no_opening_pass}; "
        f"kern_utilisation_min_max="
        f"({min(utilisations) if utilisations else 0.0:.6g},"
        f"{max(utilisations) if utilisations else 0.0:.6g}); "
        f"full_contact_pressure_min_max_MPa="
        f"({min(minimum_pressures) if minimum_pressures else 0.0:.6g},"
        f"{max(maximum_pressures) if maximum_pressures else 0.0:.6g}); "
        f"joint_area_min_max_m2="
        f"({min((record['area'] for record in records), default=0.0):.6g},"
        f"{max((record['area'] for record in records), default=0.0):.6g}); "
        f"eccentricity_min_max_m="
        f"({min((record['eccentricity'] for record in records), default=0.0):.6g},"
        f"{max((record['eccentricity'] for record in records), default=0.0):.6g}); "
        f"max_tangential_force_kN="
        f"{max((record['tangential_force'] for record in records), default=0.0):.6g}; "
        f"failures={failures}; valid={valid}; "
        "pressure_model=linear_full_contact_diagnostic; "
        "negative_min_pressure_means_opening; "
        "partial_contact_crushing_friction_and_pin_capacity_not_checked"
    )
    return (
        joint_kerns,
        utilisations,
        minimum_pressures,
        maximum_pressures,
        opening_ids,
        no_opening_pass,
        valid,
        report,
    )


def _case_key(case_index, angle, identifier):
    return "case[{0}]/angle[{1:.6g}]/{2}".format(
        case_index, angle, identifier
    )


def create_directional_joint_opening_and_pressure(
    joint_profiles,
    joint_planes,
    joint_ids,
    case_angles,
    case_joint_ids,
    contact_points,
    contact_resultants,
):
    identifiers = [str(value) for value in _as_list(joint_ids)]
    try:
        angles = [float(value) for value in _as_list(case_angles)]
    except (TypeError, ValueError):
        raise ValueError("CaseAngles must contain numeric degree values.")
    row_ids = [str(value) for value in _as_list(case_joint_ids)]
    points = _as_list(contact_points)
    resultants = _as_list(contact_resultants)
    joint_count = len(identifiers)
    expected_rows = len(angles) * joint_count
    if (
        not joint_count
        or not angles
        or len(set(identifiers)) != joint_count
        or len(row_ids) != expected_rows
        or len(points) != expected_rows
        or len(resultants) != expected_rows
    ):
        raise ValueError(
            "Stage 07C directional inputs must contain one complete "
            "JointIDs-sized block for every CaseAngles item."
        )

    profiles = [_coerce_curve(value) for value in _as_list(joint_profiles)]
    planes = [_coerce_plane(value) for value in _as_list(joint_planes)]
    if (
        len(profiles) != joint_count
        or len(planes) != joint_count
        or any(
            profile is None or not profile.IsValid or not profile.IsClosed
            for profile in profiles
        )
        or any(plane is None or not plane.IsValid for plane in planes)
    ):
        raise ValueError(
            "JointProfiles, JointPlanes and JointIDs must correspond one-to-one."
        )
    prepared_sections = []
    for profile, plane in zip(profiles, planes):
        try:
            prepared_sections.append(_prepare_joint_section(profile, plane))
        except Exception:
            prepared_sections.append(None)

    joint_kerns = []
    utilisations = []
    minimum_pressures = []
    maximum_pressures = []
    opening_keys = []
    failures = []
    for case_index, angle in enumerate(angles):
        start = case_index * joint_count
        end = start + joint_count
        if row_ids[start:end] != identifiers:
            failures.append((case_index, angle, "case_joint_id_order_mismatch"))
            continue
        result = create_joint_opening_and_pressure(
            profiles,
            planes,
            points[start:end],
            resultants[start:end],
            identifiers,
            prepared_sections,
        )
        if not result[6]:
            failures.append((case_index, angle, result[7]))
            continue
        joint_kerns.extend(result[0])
        utilisations.extend(result[1])
        minimum_pressures.extend(result[2])
        maximum_pressures.extend(result[3])
        opening_keys.extend(
            _case_key(case_index, angle, identifier)
            for identifier in result[4]
        )

    valid = (
        not failures
        and len(joint_kerns)
        == len(utilisations)
        == len(minimum_pressures)
        == len(maximum_pressures)
        == expected_rows
    )
    no_opening_pass = valid and not opening_keys
    report = (
        f"stage=directional_joint_opening_and_pressure; "
        f"angle_cases={len(angles)}; joints_per_case={joint_count}; "
        f"flattened_rows={len(utilisations)}; opening_case_rows={len(opening_keys)}; "
        f"kern_utilisation_min_max="
        f"({min(utilisations) if utilisations else 0.0:.6g},"
        f"{max(utilisations) if utilisations else 0.0:.6g}); "
        f"full_contact_pressure_min_max_MPa="
        f"({min(minimum_pressures) if minimum_pressures else 0.0:.6g},"
        f"{max(maximum_pressures) if maximum_pressures else 0.0:.6g}); "
        f"case_failures={failures}; no_opening_pass={no_opening_pass}; "
        f"valid={valid}; flattening=case_major; "
        "negative_min_pressure_means_opening"
    )
    return (
        joint_kerns,
        utilisations,
        minimum_pressures,
        maximum_pressures,
        opening_keys,
        no_opening_pass,
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
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        JointKerns,
        KernUtilisations,
        MinPressures,
        MaxPressures,
        OpeningCaseKeys,
        NoOpeningPass,
        Valid,
        Report,
    ) = create_directional_joint_opening_and_pressure(*_component_inputs)
