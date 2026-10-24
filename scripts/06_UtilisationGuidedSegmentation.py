#! python 3
"""Stage 06 - utilisation-guided radially graded stone segmentation.

Grasshopper Python inputs:
    WholeSolid         Brep         WholeSolid from Stage 01 (Item)
    TNAEdges           list[Line]   TNAEdges from Stage 04 (List)
    Utilisations       list[float]  Utilisations from Stage 05 (List)
    TargetBlockLength  float        Desired block length (Item)
    CrownBlockRadius   float        Minimum crown-block plan radius (Item)
    TopSlopePercent    float        Optional Stage 01 top slope (Item)

Unconnected numeric inputs use their defaults.

Outputs:
    Blocks, JointPlanes, JointProfiles, BlockIDs, Valid, Report

The radial TNA graph is reconstructed from TNAEdges. TargetBlockLength controls
only the approximate number of blocks. Integrated section volume is then
distributed with increasing crown-to-support weights, so blocks near the
centre are generally smaller and blocks toward the supports are generally
larger. Internal stations may shift locally toward lower Stage-05 utilisation.
Every cutting plane is normal to the local thrust direction. One connected
central crown block is retained; all remaining pieces are ordered by arch and
crown-to-support position. At zero or empty TopSlopePercent the historical
CrownBlockRadius branch-length rule is retained. At a positive slope, the
first TNA-normal joint is moved inward until its trace on the highest horizontal
top plane reaches, or conservatively stays just inside, that plan radius. This
keeps the whole centre block under the Stage 01 horizontal cap while retaining
a local discrete-TNA-normal joint; no downstream output is fed back to Stage 01.

This stage creates geometry only. Exact block weights, contact resultants,
joint opening, partial-contact pressure and sliding deficit are evaluated after
segmentation. No connector sizing or engineering strength verdict is reported.
"""

import math

import Rhino
import Rhino.Geometry as rg


DEFAULT_TARGET_BLOCK_LENGTH = 1.40
DEFAULT_CROWN_BLOCK_RADIUS = 1.60
DEFAULT_TOP_SLOPE_PERCENT = 0.0
MAX_TOP_SLOPE_PERCENT = 5.0
MAX_BLOCKS_PER_ARCH = 30
STATION_SEARCH_SAMPLES = 17
DISTANCE_PENALTY = 0.20
INNER_VOLUME_WEIGHT = 0.60
OUTER_VOLUME_WEIGHT = 1.40
JOINT_PROFILE_MATCH_SAMPLES = 24


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


def _optional_positive(value, default, name):
    if value is None:
        return default
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if number <= 0.0:
        raise ValueError("{0} must be positive.".format(name))
    return number


def _point_key(point, tolerance):
    scale = 1.0 / max(tolerance, 1.0e-9)
    return (
        int(round(point.X * scale)),
        int(round(point.Y * scale)),
        int(round(point.Z * scale)),
    )


def _radial_network(lines, tolerance):
    points = []
    key_to_index = {}
    edge_count = 0
    adjacency = []

    for line in lines:
        indices = []
        for point in (line.From, line.To):
            key = _point_key(point, tolerance)
            if key not in key_to_index:
                key_to_index[key] = len(points)
                points.append(rg.Point3d(point))
                adjacency.append([])
            indices.append(key_to_index[key])
        edge_index = edge_count
        a, b = indices
        adjacency[a].append((b, edge_index))
        adjacency[b].append((a, edge_index))
        edge_count += 1

    crown_candidates = [
        index for index, neighbours in enumerate(adjacency)
        if len(neighbours) > 2
    ]
    if len(crown_candidates) != 1:
        raise ValueError(
            "TNAEdges must contain one radial crown node."
        )
    crown = crown_candidates[0]
    branches = []
    branch_edges = []
    for neighbour, first_edge in adjacency[crown]:
        nodes = [crown, neighbour]
        edges = [first_edge]
        previous = crown
        current = neighbour
        visited = {first_edge}
        while len(adjacency[current]) == 2:
            next_items = [
                item for item in adjacency[current]
                if item[0] != previous
                and item[1] not in visited
            ]
            if not next_items:
                break
            next_node, next_edge = next_items[0]
            nodes.append(next_node)
            edges.append(next_edge)
            visited.add(next_edge)
            previous, current = current, next_node
        if len(adjacency[current]) != 1:
            raise ValueError(
                "Every radial branch must terminate at one leaf support."
            )
        branches.append(nodes)
        branch_edges.append(edges)

    if sum(len(edges) for edges in branch_edges) != edge_count:
        raise ValueError("Some TNAEdges do not belong to a radial branch.")
    return points, crown, branches, branch_edges


def _branch_data(points, branch_nodes, branch_edge_indices,
                 utilisations):
    branch_points = [points[index] for index in branch_nodes]
    lengths = []
    cumulative = [0.0]
    for start, end in zip(branch_points[:-1], branch_points[1:]):
        length = start.DistanceTo(end)
        if length <= 1.0e-9:
            raise ValueError("A branch contains a zero-length edge.")
        lengths.append(length)
        cumulative.append(cumulative[-1] + length)
    branch_utilisations = [
        utilisations[index] for index in branch_edge_indices
    ]
    return (
        branch_points,
        lengths,
        cumulative,
        branch_utilisations,
    )


def _edge_at_distance(cumulative, distance):
    if distance <= cumulative[0]:
        return 0
    for index in range(len(cumulative) - 1):
        if distance <= cumulative[index + 1] + 1.0e-12:
            return index
    return len(cumulative) - 2


def _point_tangent_at_distance(branch_points, lengths, cumulative,
                               distance):
    total = cumulative[-1]
    distance = min(max(distance, 0.0), total)
    edge_index = _edge_at_distance(cumulative, distance)
    local = (
        (distance - cumulative[edge_index]) / lengths[edge_index]
    )
    start = branch_points[edge_index]
    end = branch_points[edge_index + 1]
    point = rg.Point3d(
        start.X + (end.X - start.X) * local,
        start.Y + (end.Y - start.Y) * local,
        start.Z + (end.Z - start.Z) * local,
    )
    tangent = rg.Vector3d(
        end.X - start.X,
        end.Y - start.Y,
        end.Z - start.Z,
    )
    if not tangent.Unitize():
        raise RuntimeError("Could not evaluate a branch tangent.")
    return point, tangent


def _utilisation_at_distance(cumulative, branch_utilisations, distance):
    edge_index = _edge_at_distance(cumulative, distance)
    return branch_utilisations[edge_index]


def _desired_stations(total_length, target_length, first_station,
                      apply_target_floor):
    minimum_tail = max(0.35 * target_length, 1.0e-3)
    first = (
        max(first_station, 0.45 * target_length)
        if apply_target_floor
        else first_station
    )
    if first <= 0.0 or first >= total_length - minimum_tail:
        raise ValueError(
            "CrownBlockRadius leaves no length for outer stone blocks."
        )
    outer_length = total_length - first
    block_count = max(
        1,
        int(math.floor(outer_length / target_length + 0.5)),
    )
    block_count = min(block_count, MAX_BLOCKS_PER_ARCH)
    return first, block_count


def _optimised_stations(desired, total_length, target_length,
                        cumulative, branch_utilisations,
                        minimum_first_station):
    selected = []
    for index, target in enumerate(desired):
        if index == 0:
            selected.append(target)
            continue
        previous_target = desired[index - 1]
        next_target = (
            total_length
            if index == len(desired) - 1
            else desired[index + 1]
        )
        local_gap = min(
            target - previous_target,
            next_target - target,
        )
        search_radius = max(0.0, 0.22 * local_gap)
        lower = max(
            target - search_radius,
            minimum_first_station,
        )
        upper = min(target + search_radius, total_length)
        best = None
        for sample_index in range(STATION_SEARCH_SAMPLES):
            fraction = (
                float(sample_index)
                / float(STATION_SEARCH_SAMPLES - 1)
            )
            distance = lower + (upper - lower) * fraction
            utilisation = _utilisation_at_distance(
                cumulative, branch_utilisations, distance
            )
            score = (
                utilisation
                + DISTANCE_PENALTY
                * abs(distance - target)
                / max(target_length, 1.0e-9)
            )
            candidate = (score, distance)
            if best is None or candidate < best:
                best = candidate
        selected.append(best[1])
    return selected


def _section_plane(origin, tangent, crown, support):
    radial = rg.Vector3d(
        support.X - crown.X,
        support.Y - crown.Y,
        0.0,
    )
    if not radial.Unitize():
        raise RuntimeError("A branch has no plan direction.")
    transverse = rg.Vector3d(-radial.Y, radial.X, 0.0)
    normal = rg.Vector3d(tangent)
    normal -= transverse * rg.Vector3d.Multiply(normal, transverse)
    if not normal.Unitize():
        raise RuntimeError("A joint normal could not be evaluated.")
    in_plane = rg.Vector3d.CrossProduct(normal, transverse)
    if not in_plane.Unitize():
        raise RuntimeError("A joint plane could not be constructed.")
    plane = rg.Plane(origin, transverse, in_plane)
    if not plane.IsValid:
        raise RuntimeError("A generated joint plane is invalid.")
    angle = rg.Vector3d.VectorAngle(plane.ZAxis, tangent)
    angle = min(angle, abs(math.pi - angle))
    return plane, math.degrees(angle)


def _plan_radial_coordinate(point, crown, support):
    radial = rg.Vector3d(
        support.X - crown.X,
        support.Y - crown.Y,
        0.0,
    )
    if not radial.Unitize():
        raise RuntimeError("A branch has no plan direction.")
    offset = rg.Vector3d(
        point.X - crown.X,
        point.Y - crown.Y,
        0.0,
    )
    return rg.Vector3d.Multiply(offset, radial)


def _joint_top_trace_radius(point, tangent, crown, support,
                            top_elevation):
    plane, _ = _section_plane(point, tangent, crown, support)
    radial = rg.Vector3d(
        support.X - crown.X,
        support.Y - crown.Y,
        0.0,
    )
    if not radial.Unitize():
        raise RuntimeError("A branch has no plan direction.")
    radial_normal = rg.Vector3d.Multiply(plane.Normal, radial)
    if abs(radial_normal) <= 1.0e-12:
        raise RuntimeError(
            "A TNA-normal joint does not intersect the crown top in a "
            "stable radial direction."
        )
    axis_radius = _plan_radial_coordinate(point, crown, support)
    return (
        axis_radius
        - plane.Normal.Z
        * (top_elevation - point.Z)
        / radial_normal
    )


def _first_station_for_top_trace(
    branch_points,
    lengths,
    cumulative,
    crown_radius,
    top_elevation,
    tolerance,
):
    crown = branch_points[0]
    support = branch_points[-1]
    radial_support = _plan_radial_coordinate(support, crown, support)
    if crown_radius >= radial_support - tolerance:
        raise ValueError(
            "CrownBlockRadius must be smaller than every support plan radius."
        )

    root_tolerance = max(tolerance * 0.25, 1.0e-8)
    edge_nudge = max(tolerance * 0.05, 1.0e-8)
    best_inside = None
    for edge_index, edge_length in enumerate(lengths):
        lower = cumulative[edge_index]
        upper = cumulative[edge_index + 1]
        if edge_index > 0:
            lower = min(upper, lower + edge_nudge)
        lower_point, lower_tangent = _point_tangent_at_distance(
            branch_points, lengths, cumulative, lower
        )
        upper_point, upper_tangent = _point_tangent_at_distance(
            branch_points, lengths, cumulative, upper
        )
        lower_value = _joint_top_trace_radius(
            lower_point,
            lower_tangent,
            crown,
            support,
            top_elevation,
        ) - crown_radius
        upper_value = _joint_top_trace_radius(
            upper_point,
            upper_tangent,
            crown,
            support,
            top_elevation,
        ) - crown_radius

        for candidate_distance, candidate_value in (
            (lower, lower_value),
            (upper, upper_value),
        ):
            if candidate_value <= root_tolerance:
                candidate = (
                    abs(candidate_value),
                    -candidate_distance,
                    candidate_distance,
                )
                if best_inside is None or candidate < best_inside:
                    best_inside = candidate

        if abs(lower_value) <= root_tolerance:
            station = lower
        elif abs(upper_value) <= root_tolerance:
            station = upper
        elif lower_value * upper_value > 0.0:
            continue
        else:
            left = lower
            right = upper
            left_value = lower_value
            for _iteration in range(64):
                middle = 0.5 * (left + right)
                point, tangent = _point_tangent_at_distance(
                    branch_points, lengths, cumulative, middle
                )
                middle_value = _joint_top_trace_radius(
                    point,
                    tangent,
                    crown,
                    support,
                    top_elevation,
                ) - crown_radius
                if abs(middle_value) <= root_tolerance:
                    left = right = middle
                    break
                if left_value * middle_value <= 0.0:
                    right = middle
                else:
                    left = middle
                    left_value = middle_value
            station = 0.5 * (left + right)

        if station <= max(tolerance * 10.0, 1.0e-6):
            raise ValueError(
                "CrownBlockRadius leaves no positive TNA-normal crown joint."
            )
        point, tangent = _point_tangent_at_distance(
            branch_points, lengths, cumulative, station
        )
        achieved = _joint_top_trace_radius(
            point,
            tangent,
            crown,
            support,
            top_elevation,
        )
        if abs(achieved - crown_radius) <= max(tolerance, 1.0e-7):
            return station

    if best_inside is not None:
        station = best_inside[2]
        if station > max(tolerance * 10.0, 1.0e-6):
            return station

    raise ValueError(
        "No TNA-normal first joint can keep its crown-top trace at or "
        "inside CrownBlockRadius for this branch."
    )


def _section_curves(solid, plane, tolerance):
    result = rg.Intersect.Intersection.BrepPlane(
        solid, plane, tolerance
    )
    if not result or not result[0]:
        return []
    raw_curves = list(result[1]) if result[1] else []
    return [
        curve for curve in rg.Curve.JoinCurves(raw_curves, tolerance)
        if curve is not None and curve.IsValid and curve.IsClosed
    ]


def _curve_area(curve):
    properties = rg.AreaMassProperties.Compute(curve)
    if properties is None:
        return 0.0
    area = abs(properties.Area)
    properties.Dispose()
    return area


def _containing_profile(curves, point, plane, tolerance):
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
            "No joint profile containing the TNA station was found."
        )
    profile = max(containing, key=_curve_area)
    for curve in curves:
        if curve is not profile:
            curve.Dispose()
    return profile


def _distance_at_accumulated_value(distances, accumulated, target):
    if target <= accumulated[0]:
        return distances[0]
    for index in range(len(accumulated) - 1):
        upper = accumulated[index + 1]
        if target <= upper + 1.0e-12:
            lower = accumulated[index]
            interval = upper - lower
            if interval <= 1.0e-12:
                return distances[index]
            fraction = (target - lower) / interval
            return (
                distances[index]
                + (distances[index + 1] - distances[index])
                * fraction
            )
    return distances[-1]


def _graded_volume_stations(
    solid,
    branch_points,
    lengths,
    cumulative,
    target_length,
    crown_radius,
    top_elevation,
    tolerance,
    align_plan_radius,
):
    crown = branch_points[0]
    plan_radii = [
        math.hypot(point.X - crown.X, point.Y - crown.Y)
        for point in branch_points
    ]
    if any(
        plan_radii[index + 1] <= plan_radii[index] + tolerance
        for index in range(len(plan_radii) - 1)
    ):
        raise ValueError(
            "A branch must move monotonically outward in plan."
        )
    if align_plan_radius:
        first_station = _first_station_for_top_trace(
            branch_points,
            lengths,
            cumulative,
            crown_radius,
            top_elevation,
            tolerance,
        )
    else:
        first_station = crown_radius
    first, block_count = _desired_stations(
        cumulative[-1], target_length, first_station,
        not align_plan_radius,
    )
    outer_length = cumulative[-1] - first
    evaluation_end = cumulative[-1] - max(
        tolerance * 10.0,
        outer_length * 1.0e-5,
    )
    if evaluation_end <= first:
        raise ValueError(
            "The outer branch is too short for volume integration."
        )
    sample_count = max(49, min(161, 8 * len(lengths) + 1))
    distances = [
        first
        + (evaluation_end - first)
        * index
        / float(sample_count - 1)
        for index in range(sample_count)
    ]
    areas = []
    for distance in distances:
        point, tangent = _point_tangent_at_distance(
            branch_points, lengths, cumulative, distance
        )
        plane, _ = _section_plane(
            point,
            tangent,
            branch_points[0],
            branch_points[-1],
        )
        profile = _containing_profile(
            _section_curves(solid, plane, tolerance),
            point,
            plane,
            tolerance,
        )
        areas.append(_curve_area(profile))
        profile.Dispose()

    accumulated = [0.0]
    for index in range(len(distances) - 1):
        step = distances[index + 1] - distances[index]
        increment = 0.5 * (areas[index] + areas[index + 1]) * step
        accumulated.append(accumulated[-1] + increment)
    estimated_volume = accumulated[-1]
    if block_count == 1:
        volume_weights = [1.0]
    else:
        volume_weights = [
            INNER_VOLUME_WEIGHT
            + (OUTER_VOLUME_WEIGHT - INNER_VOLUME_WEIGHT)
            * block_index
            / float(block_count - 1)
            for block_index in range(block_count)
        ]
    weight_total = sum(volume_weights)
    stations = [first]
    accumulated_weight = 0.0
    for block_index in range(1, block_count):
        accumulated_weight += volume_weights[block_index - 1]
        target_volume = (
            estimated_volume * accumulated_weight / weight_total
        )
        stations.append(
            _distance_at_accumulated_value(
                distances, accumulated, target_volume
            )
        )
    return stations, block_count


def _plane_coordinates(point, plane):
    vector = rg.Vector3d(
        point.X - plane.Origin.X,
        point.Y - plane.Origin.Y,
        point.Z - plane.Origin.Z,
    )
    return (
        rg.Vector3d.Multiply(vector, plane.XAxis),
        rg.Vector3d.Multiply(vector, plane.YAxis),
    )


def _local_cutter(profile, plane, tolerance):
    parameters = profile.DivideByCount(96, True)
    if not parameters:
        raise RuntimeError("A joint profile could not be sampled.")
    coordinates = [
        _plane_coordinates(profile.PointAt(parameter), plane)
        for parameter in parameters
    ]
    minimum_x = min(point[0] for point in coordinates)
    maximum_x = max(point[0] for point in coordinates)
    minimum_y = min(point[1] for point in coordinates)
    maximum_y = max(point[1] for point in coordinates)
    size = max(
        maximum_x - minimum_x,
        maximum_y - minimum_y,
        tolerance * 10.0,
    )
    margin = max(size * 0.25, tolerance * 20.0)
    surface = rg.PlaneSurface(
        plane,
        rg.Interval(minimum_x - margin, maximum_x + margin),
        rg.Interval(minimum_y - margin, maximum_y + margin),
    )
    cutter = surface.ToBrep()
    surface.Dispose()
    if cutter is None or not cutter.IsValid:
        raise RuntimeError("A local joint cutter is invalid.")
    return cutter


def _normalise_solid(brep):
    if brep.SolidOrientation == rg.BrepSolidOrientation.Inward:
        brep.Flip()
    return brep.IsValid and brep.IsSolid


def _close_split_piece(brep, tolerance):
    if brep is None or not brep.IsValid:
        return None
    candidate = brep
    if not candidate.IsSolid:
        capped = candidate.CapPlanarHoles(tolerance)
        if capped is not None:
            candidate = capped
    return candidate if _normalise_solid(candidate) else None


def _split_all(solid, cutters, tolerance):
    pieces = [solid.DuplicateBrep()]
    for cutter_index, cutter in enumerate(cutters):
        next_pieces = []
        split_count = 0
        for piece in pieces:
            result = piece.Split(cutter, tolerance)
            candidates = list(result) if result is not None else []
            valid = []
            for candidate in candidates:
                closed = _close_split_piece(candidate, tolerance)
                if closed is not None:
                    valid.append(closed)
            if len(valid) >= 2:
                next_pieces.extend(valid)
                split_count += 1
            else:
                next_pieces.append(piece)
        if split_count != 1:
            raise RuntimeError(
                "Joint cutter {0} split {1} pieces; expected exactly one."
                .format(cutter_index, split_count)
            )
        pieces = next_pieces
    return pieces


def _volume_centroid(solid):
    properties = rg.VolumeMassProperties.Compute(solid)
    if properties is None:
        raise RuntimeError("Could not evaluate a block volume.")
    volume = abs(properties.Volume)
    centroid = properties.Centroid
    properties.Dispose()
    return volume, centroid


def _connected_component_count(brep):
    result = brep.GetConnectedComponents()
    components = list(result) if result else []
    if not components:
        return 1
    count = len(components)
    for component in components:
        component.Dispose()
    return count


def _area_centroid(geometry, description):
    properties = rg.AreaMassProperties.Compute(geometry)
    if properties is None:
        raise RuntimeError(
            "Could not evaluate {0} area properties.".format(description)
        )
    area = abs(properties.Area)
    centroid = rg.Point3d(properties.Centroid)
    properties.Dispose()
    if area <= 0.0 or not centroid.IsValid:
        raise RuntimeError(
            "{0} has invalid area properties.".format(description)
        )
    return area, centroid


def _center_top_flatness(block, source_solid, top_elevation, tolerance):
    upward_face_count = 0
    minimum_top_vertex_z = top_elevation
    maximum_top_drop = 0.0
    face_tolerance = max(5.0 * tolerance, 1.0e-7)
    minimum_top_normal_z = (
        1.0 / math.sqrt(
            1.0 + (MAX_TOP_SLOPE_PERCENT / 100.0) ** 2
        )
        - 1.0e-6
    )
    boundary_probe = max(10.0 * tolerance, 1.0e-4)
    for face in block.Faces:
        domain_u = face.Domain(0)
        domain_v = face.Domain(1)
        normal = face.NormalAt(
            0.5 * (domain_u.T0 + domain_u.T1),
            0.5 * (domain_v.T0 + domain_v.T1),
        )
        if not normal.Unitize() or abs(normal.Z) < minimum_top_normal_z:
            continue
        face_brep = face.DuplicateFace(False)
        if face_brep is None or not face_brep.IsValid:
            if face_brep is not None:
                face_brep.Dispose()
            raise RuntimeError(
                "Could not inspect an upward-facing crown-block face."
            )
        try:
            elevations = [
                vertex.Location.Z for vertex in face_brep.Vertices
            ]
            _face_area, face_centroid = _area_centroid(
                face_brep, "an upward crown-block face"
            )
        finally:
            face_brep.Dispose()
        if not elevations:
            raise RuntimeError(
                "An upward-facing crown-block face has no vertices."
            )
        if max(elevations) < top_elevation - face_tolerance:
            continue
        plus_point = rg.Point3d(
            face_centroid.X + normal.X * boundary_probe,
            face_centroid.Y + normal.Y * boundary_probe,
            face_centroid.Z + normal.Z * boundary_probe,
        )
        minus_point = rg.Point3d(
            face_centroid.X - normal.X * boundary_probe,
            face_centroid.Y - normal.Y * boundary_probe,
            face_centroid.Z - normal.Z * boundary_probe,
        )
        plus_inside = source_solid.IsPointInside(
            plus_point, tolerance, True
        )
        minus_inside = source_solid.IsPointInside(
            minus_point, tolerance, True
        )
        if plus_inside == minus_inside:
            continue
        upward_face_count += 1
        face_minimum = min(elevations)
        minimum_top_vertex_z = min(
            minimum_top_vertex_z, face_minimum
        )
        maximum_top_drop = max(
            maximum_top_drop,
            top_elevation - face_minimum,
        )
    return {
        "horizontal": (
            upward_face_count > 0
            and maximum_top_drop <= face_tolerance
        ),
        "upward_face_count": upward_face_count,
        "minimum_top_vertex_z": minimum_top_vertex_z,
        "maximum_top_drop": maximum_top_drop,
        "tolerance": face_tolerance,
    }


def _planar_block_faces(blocks, tolerance):
    records = []
    for block_index, block in enumerate(blocks):
        for face_index, face in enumerate(block.Faces):
            if not face.IsPlanar(tolerance):
                continue
            face_brep = face.DuplicateFace(False)
            if face_brep is None or not face_brep.IsValid:
                if face_brep is not None:
                    face_brep.Dispose()
                continue
            try:
                area, centroid = _area_centroid(
                    face_brep, "a planar block face"
                )
                domain_u = face.Domain(0)
                domain_v = face.Domain(1)
                normal = face.NormalAt(
                    0.5 * (domain_u.T0 + domain_u.T1),
                    0.5 * (domain_v.T0 + domain_v.T1),
                )
                if not normal.Unitize():
                    raise RuntimeError(
                        "A planar block face has no valid normal."
                    )
            except Exception:
                face_brep.Dispose()
                raise
            records.append(
                (
                    block_index,
                    face_index,
                    face_brep,
                    area,
                    centroid,
                    normal,
                )
            )
    return records


def _profile_sample_points(profile):
    parameters = profile.DivideByCount(
        JOINT_PROFILE_MATCH_SAMPLES, True
    )
    if not parameters:
        raise RuntimeError("A joint profile could not be sampled.")
    return [profile.PointAt(parameter) for parameter in parameters]


def _joint_profile_data(profile, tolerance):
    profile_area, profile_centroid = _area_centroid(
        profile, "a joint profile"
    )
    characteristic_length = max(math.sqrt(profile_area), tolerance)
    return (
        profile_area,
        profile_centroid,
        _profile_sample_points(profile),
        max(5.0 * tolerance, characteristic_length * 1.0e-6),
        max(5.0 * tolerance, characteristic_length * 1.0e-5),
        max(5.0 * tolerance, characteristic_length * 1.0e-5),
        max(
            profile_area * 1.0e-4,
            25.0 * tolerance * tolerance,
        ),
    )


def _profile_face_match(profile_data, plane, face_record):
    (
        block_index,
        face_index,
        face_brep,
        face_area,
        face_centroid,
        face_normal,
    ) = face_record
    (
        profile_area,
        profile_centroid,
        sample_points,
        plane_tolerance,
        centroid_tolerance,
        boundary_tolerance,
        area_tolerance,
    ) = profile_data
    normal_alignment = abs(
        rg.Vector3d.Multiply(face_normal, plane.Normal)
    )
    if normal_alignment < math.cos(math.radians(0.01)):
        return None
    plane_error = abs(plane.DistanceTo(face_centroid))
    if plane_error > plane_tolerance:
        return None
    area_error = abs(face_area - profile_area)
    if area_error > area_tolerance:
        return None
    centroid_error = face_centroid.DistanceTo(profile_centroid)
    if centroid_error > centroid_tolerance:
        return None
    boundary_error = 0.0
    for point in sample_points:
        closest = face_brep.ClosestPoint(point)
        if closest is None or not closest.IsValid:
            return None
        boundary_error = max(boundary_error, point.DistanceTo(closest))
    if boundary_error > boundary_tolerance:
        return None
    return (
        block_index,
        face_index,
        area_error,
        centroid_error,
        boundary_error,
        plane_error,
        face_normal,
    )


def _validate_joint_profiles_against_blocks(
    blocks, joint_planes, joint_profiles, tolerance
):
    face_records = _planar_block_faces(blocks, tolerance)
    failures = []
    used_faces = set()
    reused_faces = []
    maximum_area_error = 0.0
    maximum_centroid_error = 0.0
    maximum_boundary_error = 0.0
    maximum_plane_error = 0.0
    matched_profiles = 0
    try:
        for joint_index, (plane, profile) in enumerate(
            zip(joint_planes, joint_profiles)
        ):
            profile_data = _joint_profile_data(profile, tolerance)
            matches = []
            for record in face_records:
                match = _profile_face_match(
                    profile_data, plane, record
                )
                if match is not None:
                    matches.append(match)
            matched_blocks = sorted(
                set(match[0] for match in matches)
            )
            opposite_normals = False
            if len(matches) == 2 and len(matched_blocks) == 2:
                opposite_normals = (
                    rg.Vector3d.Multiply(
                        matches[0][6], matches[1][6]
                    )
                    <= -math.cos(math.radians(0.01))
                )
            if (
                len(matches) != 2
                or len(matched_blocks) != 2
                or not opposite_normals
            ):
                failures.append(
                    (
                        joint_index,
                        len(matches),
                        matched_blocks,
                        opposite_normals,
                    )
                )
                continue
            face_keys = [
                (match[0], match[1]) for match in matches
            ]
            repeated = [
                key for key in face_keys if key in used_faces
            ]
            if repeated:
                reused_faces.extend(
                    (joint_index, key) for key in repeated
                )
                failures.append(
                    (
                        joint_index,
                        len(matches),
                        matched_blocks,
                        "reused_face",
                    )
                )
                continue
            used_faces.update(face_keys)
            matched_profiles += 1
            maximum_area_error = max(
                maximum_area_error,
                max(match[2] for match in matches),
            )
            maximum_centroid_error = max(
                maximum_centroid_error,
                max(match[3] for match in matches),
            )
            maximum_boundary_error = max(
                maximum_boundary_error,
                max(match[4] for match in matches),
            )
            maximum_plane_error = max(
                maximum_plane_error,
                max(match[5] for match in matches),
            )
    finally:
        for record in face_records:
            record[2].Dispose()
    return {
        "valid": (
            matched_profiles == len(joint_profiles)
            and not failures
            and not reused_faces
        ),
        "matched_profiles": matched_profiles,
        "failures": failures,
        "reused_faces": reused_faces,
        "maximum_area_error": maximum_area_error,
        "maximum_centroid_error": maximum_centroid_error,
        "maximum_boundary_error": maximum_boundary_error,
        "maximum_plane_error": maximum_plane_error,
    }


def _closest_branch(point, branch_curves):
    best = None
    for branch_index, curve in enumerate(branch_curves):
        result = curve.ClosestPoint(point)
        if not result or not result[0]:
            continue
        parameter = result[1]
        closest = curve.PointAt(parameter)
        distance = closest.DistanceTo(point)
        candidate = (distance, branch_index, parameter)
        if best is None or candidate < best:
            best = candidate
    if best is None:
        raise RuntimeError("A block could not be assigned to a branch.")
    return best[1], best[2]


def _order_and_label_blocks(pieces, crown_point, branch_curves,
                            tolerance):
    records = []
    center_candidates = []
    for piece_index, piece in enumerate(pieces):
        volume, centroid = _volume_centroid(piece)
        if volume <= 0.0:
            raise RuntimeError("A split block has non-positive volume.")
        if piece.IsPointInside(crown_point, tolerance, False):
            center_candidates.append(piece_index)
        branch_index, parameter = _closest_branch(
            centroid, branch_curves
        )
        records.append(
            (piece_index, piece, centroid, branch_index, parameter)
        )
    if len(center_candidates) == 1:
        center_index = center_candidates[0]
    else:
        center_index = min(
            records,
            key=lambda item: math.hypot(
                item[2].X - crown_point.X,
                item[2].Y - crown_point.Y,
            ),
        )[0]

    ordered_blocks = []
    block_ids = []
    center_record = next(
        record for record in records if record[0] == center_index
    )
    ordered_blocks.append(center_record[1])
    block_ids.append("center/block[0]")

    for branch_index in range(len(branch_curves)):
        branch_records = [
            record for record in records
            if record[0] != center_index
            and record[3] == branch_index
        ]
        branch_records.sort(key=lambda item: item[4])
        for local_index, record in enumerate(branch_records):
            ordered_blocks.append(record[1])
            block_ids.append(
                "arch[{0}]/block[{1}]".format(
                    branch_index + 1, local_index + 1
                )
            )
    if len(ordered_blocks) != len(pieces):
        raise RuntimeError("Not every split block received a stable ID.")
    return (
        ordered_blocks,
        block_ids,
        len(center_candidates),
    )


def _radial_volume_trends(block_ids, block_volumes):
    """Summarise the overall inner-to-outer volume trend per branch."""
    grouped = {}
    for block_id, volume in zip(block_ids, block_volumes):
        if not block_id.startswith("arch["):
            continue
        arch_index = int(block_id.split("[", 1)[1].split("]", 1)[0])
        grouped.setdefault(arch_index, []).append(float(volume))

    records = []
    for arch_index in sorted(grouped):
        volumes = grouped[arch_index]
        if len(volumes) == 1:
            records.append(
                (
                    arch_index,
                    1,
                    volumes[0],
                    volumes[0],
                    1.0,
                    False,
                )
            )
            continue
        bin_count = max(1, int(math.ceil(len(volumes) / 3.0)))
        inner_mean = sum(volumes[:bin_count]) / float(bin_count)
        outer_mean = sum(volumes[-bin_count:]) / float(bin_count)
        outer_inner_ratio = outer_mean / max(inner_mean, 1.0e-12)
        records.append(
            (
                arch_index,
                len(volumes),
                inner_mean,
                outer_mean,
                outer_inner_ratio,
                outer_mean > inner_mean,
            )
        )
    return records


def create_utilisation_guided_segmentation(
    whole_solid,
    tna_edges,
    utilisations,
    target_block_length=None,
    crown_block_radius=None,
    top_slope_percent=None,
):
    solid = _coerce_brep(whole_solid)
    if solid is None or not solid.IsValid or not solid.IsSolid:
        raise ValueError("WholeSolid must be one valid closed Brep.")
    lines = []
    for edge_index, value in enumerate(_as_list(tna_edges)):
        line = _coerce_line(value)
        if line is None or not line.IsValid or line.Length <= 0.0:
            raise ValueError(
                "TNAEdges contains an invalid line at index {0}."
                .format(edge_index)
            )
        lines.append(line)
    if not lines:
        raise ValueError("TNAEdges must contain valid lines.")
    utilisation_values = [
        float(value) for value in _as_list(utilisations)
    ]
    if len(utilisation_values) != len(lines):
        raise ValueError(
            "Utilisations must correspond one-to-one with TNAEdges."
        )
    if any(
        not math.isfinite(value) or value < 0.0
        for value in utilisation_values
    ):
        raise ValueError(
            "Utilisations must contain finite non-negative numbers."
        )
    target_length = _optional_positive(
        target_block_length,
        DEFAULT_TARGET_BLOCK_LENGTH,
        "TargetBlockLength",
    )
    crown_radius = _optional_positive(
        crown_block_radius,
        DEFAULT_CROWN_BLOCK_RADIUS,
        "CrownBlockRadius",
    )
    top_slope_percent = (
        DEFAULT_TOP_SLOPE_PERCENT
        if top_slope_percent is None
        else float(top_slope_percent)
    )
    if (
        not math.isfinite(top_slope_percent)
        or top_slope_percent < 0.0
        or top_slope_percent > MAX_TOP_SLOPE_PERCENT
    ):
        raise ValueError(
            "TopSlopePercent must be finite and between 0 and 5."
        )
    align_plan_radius = top_slope_percent > 0.0
    tolerance = _document_tolerance()
    (
        points,
        crown,
        branches,
        branch_edges,
    ) = _radial_network(
        lines, max(tolerance * 0.1, 1.0e-7)
    )

    joint_planes = []
    joint_profiles = []
    cutters = []
    first_station_distances = []
    first_station_plan_radii = []
    per_arch_blocks = []
    high_risk_joints = 0
    shifted_joint_count = 0
    maximum_station_shift = 0.0
    maximum_angle_error = 0.0
    branch_curves = []
    first_station_top_trace_radii = []
    first_station_top_trace_gaps = []
    top_elevation = solid.GetBoundingBox(True).Max.Z

    for branch_nodes, edge_indices in zip(branches, branch_edges):
        (
            branch_points,
            lengths,
            cumulative,
            branch_utilisations,
        ) = _branch_data(
            points,
            branch_nodes,
            edge_indices,
            utilisation_values,
        )
        branch_curves.append(rg.PolylineCurve(branch_points))
        desired, outer_block_count = _graded_volume_stations(
            solid,
            branch_points,
            lengths,
            cumulative,
            target_length,
            crown_radius,
            top_elevation,
            tolerance,
            align_plan_radius,
        )
        selected = _optimised_stations(
            desired,
            cumulative[-1],
            target_length,
            cumulative,
            branch_utilisations,
            desired[0],
        )
        first_station_distances.append(selected[0])
        first_point, first_tangent = _point_tangent_at_distance(
            branch_points, lengths, cumulative, selected[0]
        )
        first_station_plan_radii.append(
            math.hypot(
                first_point.X - branch_points[0].X,
                first_point.Y - branch_points[0].Y,
            )
        )
        first_station_top_trace_radii.append(
            _joint_top_trace_radius(
                first_point,
                first_tangent,
                branch_points[0],
                branch_points[-1],
                top_elevation,
            )
        )
        first_station_top_trace_gaps.append(
            crown_radius - first_station_top_trace_radii[-1]
        )
        per_arch_blocks.append(outer_block_count)
        high_risk_joints += sum(
            _utilisation_at_distance(
                cumulative, branch_utilisations, distance
            ) > 1.0
            for distance in selected
        )
        for desired_distance, selected_distance in zip(
            desired, selected
        ):
            shift = abs(selected_distance - desired_distance)
            shifted_joint_count += shift > tolerance
            maximum_station_shift = max(
                maximum_station_shift, shift
            )

        for distance in selected:
            point, tangent = _point_tangent_at_distance(
                branch_points, lengths, cumulative, distance
            )
            plane, angle_error = _section_plane(
                point,
                tangent,
                branch_points[0],
                branch_points[-1],
            )
            profile = _containing_profile(
                _section_curves(solid, plane, tolerance),
                point,
                plane,
                tolerance,
            )
            cutter = _local_cutter(profile, plane, tolerance)
            joint_planes.append(plane)
            joint_profiles.append(profile)
            cutters.append(cutter)
            maximum_angle_error = max(
                maximum_angle_error, angle_error
            )

    try:
        pieces = _split_all(solid, cutters, tolerance)
    finally:
        for cutter in cutters:
            cutter.Dispose()
    expected_blocks = 1 + sum(per_arch_blocks)
    if len(pieces) != expected_blocks:
        raise RuntimeError(
            "Segmentation produced {0} blocks; expected {1}."
            .format(len(pieces), expected_blocks)
        )
    crown_point = points[crown]
    try:
        (
            blocks,
            block_ids,
            center_candidate_count,
        ) = _order_and_label_blocks(
            pieces, crown_point, branch_curves, tolerance
        )
    finally:
        for branch_curve in branch_curves:
            branch_curve.Dispose()
    center_selection_exact = center_candidate_count == 1
    center_top_flatness = (
        _center_top_flatness(
            blocks[0], solid, top_elevation, tolerance
        )
        if align_plan_radius
        else {
            "horizontal": True,
            "upward_face_count": 0,
            "minimum_top_vertex_z": top_elevation,
            "maximum_top_drop": 0.0,
            "tolerance": max(5.0 * tolerance, 1.0e-7),
        }
    )
    if align_plan_radius and not center_top_flatness["horizontal"]:
        raise RuntimeError(
            "The centre crown block exposes a sloped top face; no valid "
            "horizontal crown top was produced: {0}."
            .format(center_top_flatness)
        )

    source_volume, _ = _volume_centroid(solid)
    block_volumes = [
        _volume_centroid(block)[0] for block in blocks
    ]
    split_volume = sum(block_volumes)
    volume_error = abs(split_volume - source_volume)
    allowed_volume_error = max(
        source_volume * 2.0e-5,
        100.0 * tolerance ** 3,
    )
    radial_volume_trends = _radial_volume_trends(
        block_ids, block_volumes
    )
    radial_gradient_pass = (
        len(radial_volume_trends) == len(branches)
        and all(
            record[0] == branch_index + 1
            and record[1] == per_arch_blocks[branch_index]
            and record[5]
            for branch_index, record in enumerate(radial_volume_trends)
        )
    )
    disconnected_block_ids = [
        block_id
        for block_id, block in zip(block_ids, blocks)
        if _connected_component_count(block) != 1
    ]
    block_ids_unique = len(set(block_ids)) == len(block_ids)
    joint_face_validation = _validate_joint_profiles_against_blocks(
        blocks, joint_planes, joint_profiles, tolerance
    )
    valid = (
        len(blocks) == expected_blocks
        and len(block_ids) == expected_blocks
        and len(joint_planes) == sum(per_arch_blocks)
        and len(joint_profiles) == len(joint_planes)
        and all(block.IsValid and block.IsSolid for block in blocks)
        and all(
            profile.IsValid and profile.IsClosed
            for profile in joint_profiles
        )
        and maximum_angle_error <= 1.0e-5
        and volume_error <= allowed_volume_error
        and (
            all(
                radius <= crown_radius + tolerance
                for radius in first_station_top_trace_radii
            )
            if align_plan_radius
            else all(
                distance >= crown_radius - tolerance
                for distance in first_station_distances
            )
        )
        and radial_gradient_pass
        and center_selection_exact
        and center_top_flatness["horizontal"]
        and not disconnected_block_ids
        and block_ids_unique
        and joint_face_validation["valid"]
    )
    reported_trends = [
        (
            record[0],
            record[1],
            round(record[2], 6),
            round(record[3], 6),
            round(record[4], 6),
            record[5],
        )
        for record in radial_volume_trends
    ]
    report = (
        f"stage=utilisation_guided_segmentation; arches={len(branches)}; "
        f"target_block_length={target_length:.6g}; "
        f"crown_block_radius={crown_radius:.6g}; "
        f"top_slope_percent={top_slope_percent:.6g}; "
        "crown_radius_mode="
        f"{'joint_top_trace_at_or_inside_slope_break' if align_plan_radius else 'branch_length_legacy'}; "
        f"outer_blocks_per_arch={per_arch_blocks}; "
        f"center_block_output=1; blocks={len(blocks)}; "
        f"joints={len(joint_planes)}; "
        f"high_risk_joints={high_risk_joints}; "
        f"shifted_joints={shifted_joint_count}; "
        f"max_station_shift={maximum_station_shift:.6g}; "
        "radial_volume_trends="
        "(arch,count,inner_mean,outer_mean,outer_inner_ratio,"
        f"pass){reported_trends}; "
        f"radial_gradient_pass={radial_gradient_pass}; "
        "first_station_distances="
        f"{[round(value, 6) for value in first_station_distances]}; "
        "first_station_plan_radii="
        f"{[round(value, 6) for value in first_station_plan_radii]}; "
        "first_station_top_trace_radii="
        f"{[round(value, 6) for value in first_station_top_trace_radii]}; "
        "first_station_top_trace_gaps="
        f"{[round(value, 6) for value in first_station_top_trace_gaps]}; "
        f"center_top_horizontal={center_top_flatness['horizontal']}; "
        "center_top_upward_faces="
        f"{center_top_flatness['upward_face_count']}; "
        "center_top_maximum_drop="
        f"{center_top_flatness['maximum_top_drop']:.6g}; "
        f"max_normal_alignment_error_deg={maximum_angle_error:.6g}; "
        f"source_volume={source_volume:.6g}; "
        f"split_volume={split_volume:.6g}; "
        f"volume_error={volume_error:.6g}; tolerance={tolerance:.6g}; "
        f"center_candidate_count={center_candidate_count}; "
        f"center_selection_exact={center_selection_exact}; "
        f"disconnected_block_ids={disconnected_block_ids}; "
        f"block_ids_unique={block_ids_unique}; "
        "joint_profiles_matched_to_two_final_faces="
        f"{joint_face_validation['matched_profiles']}/"
        f"{len(joint_profiles)}; "
        "joint_profile_match_failures="
        f"{joint_face_validation['failures']}; "
        f"reused_joint_faces={joint_face_validation['reused_faces']}; "
        "max_joint_face_area_error="
        f"{joint_face_validation['maximum_area_error']:.6g}; "
        "max_joint_face_centroid_error="
        f"{joint_face_validation['maximum_centroid_error']:.6g}; "
        "max_joint_face_boundary_error="
        f"{joint_face_validation['maximum_boundary_error']:.6g}; "
        "max_joint_face_plane_error="
        f"{joint_face_validation['maximum_plane_error']:.6g}; "
        f"valid={valid}"
    )
    return (
        blocks,
        joint_planes,
        joint_profiles,
        block_ids,
        valid,
        report,
    )


try:
    _top_slope_component_input = TopSlopePercent
except NameError:
    _top_slope_component_input = None

try:
    _component_inputs = (
        WholeSolid,
        TNAEdges,
        Utilisations,
        TargetBlockLength,
        CrownBlockRadius,
        _top_slope_component_input,
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        Blocks,
        JointPlanes,
        JointProfiles,
        BlockIDs,
        Valid,
        Report,
    ) = create_utilisation_guided_segmentation(*_component_inputs)
