#! python 3
"""Stage 07B - unified ideal and horizontal-study joint resultants.

Grasshopper Python inputs:
    TNAEdges             list[Line]   TNAEdges from Stage 04 (List)
    NodalLoads           list[Vector] NodalLoads from Stage 02 (List, kN)
    JointPlanes          list[Plane]  JointPlanes from Stage 06 (List)
    JointProfiles        list[Curve]  JointProfiles from Stage 06 (List)
    JointIDs             list[str]    JointIDs from Stage 07A (List)
    JointBlockPairs      list[pair]   JointBlockPairs from Stage 07A (List)
    HorizontalLoadRatio  float        k_h = a_h/g (Item, non-negative)
    HorizontalAngleStepDegrees float  Coarse sweep step (Item, degrees)
    RefinementAngleStepDegrees float  Local refinement step (Item, degrees)
    RefinementHalfWidthDegrees float  Initial half-width (Item, degrees)
    FrictionCoefficient float         Same empirical value as Stages 07E/07F
    SlidingSafetyFactor float         Same empirical value as Stages 07E/07F

Outputs:
    CaseAngles, CaseJointIDs, ContactPoints, ContactResultants,
    TangentialForces, Valid, Report

HorizontalLoadRatio is the dimensionless coefficient k_h = a_h/g. CaseAngles
denote the direction of the equivalent horizontal inertia force applied to the
Stage 02 lumped nodal masses; they do not denote the direction of base/ground
acceleration. The corresponding base-acceleration direction is
(CaseAngle + 180 degrees) modulo 360 degrees.

HorizontalLoadRatio=0 uses one exact axial TNA baseline case and returns
strictly zero tangential force. A positive ratio first sweeps 0 <= angle < 360
at HorizontalAngleStepDegrees. For every joint it then refines the coarse
local-maximum candidates of the signed Stage 07F sliding deficit within
plus/minus RefinementHalfWidthDegrees at RefinementAngleStepDegrees. Refining
all coarse local maxima avoids assuming that the response has only one smooth
peak. If a maximum remains on a refinement boundary, that boundary is extended
until the maximum is interior or the complete 360-degree domain has been
covered. Outputs are flattened case-major lists: each CaseAngles entry owns one
consecutive JointIDs-sized block in CaseJointIDs and the three result lists.

Only Stage 02 weights at free TNA nodes enter the internal vector-equilibrium
target. Stage 02 weights located at support nodes are reported separately and
are interpreted as balanced by local support reactions; they are not hidden or
claimed to enter the free-node solve. This model enforces translational force
equilibrium at free TNA nodes only. It does not solve individual block moments,
support-contact redistribution or dynamic response.

The method is a mass-proportional, full-direction equivalent-static inertia
sensitivity study under uniform base acceleration. It is not an earthquake
safety check, response-spectrum analysis, time-history analysis or structural
safety proof.
"""

import math

import numpy as np
import Rhino
import Rhino.Geometry as rg


COMPRESSION_EPSILON = 1.0e-8
MAX_ALIGNMENT_ERROR_DEGREES = 0.001
REFINEMENT_SCORE_TOLERANCE_KN = 1.0e-10
LOAD_COMPONENT_TOLERANCE_KN = 1.0e-9
WEIGHT_PARTITION_RELATIVE_TOLERANCE = 1.0e-12


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


def _coerce_line(value):
    if isinstance(value, rg.Line):
        return rg.Line(value.From, value.To)
    if isinstance(value, rg.LineCurve):
        return rg.Line(value.PointAtStart, value.PointAtEnd)
    if hasattr(value, "Value"):
        return _coerce_line(value.Value)
    return None


def _coerce_vector(value):
    if isinstance(value, rg.Vector3d):
        return rg.Vector3d(value)
    if hasattr(value, "Value"):
        return _coerce_vector(value.Value)
    try:
        return rg.Vector3d(value)
    except Exception:
        return None


def _validated_stage02_weights(loads):
    weights = []
    for index, load in enumerate(loads):
        components = (float(load.X), float(load.Y), float(load.Z))
        if not all(math.isfinite(value) for value in components):
            raise ValueError(
                "NodalLoads contains a non-finite vector at index {0}."
                .format(index)
            )
        component_tolerance = max(
            LOAD_COMPONENT_TOLERANCE_KN,
            abs(components[2]) * WEIGHT_PARTITION_RELATIVE_TOLERANCE,
        )
        if (
            abs(components[0]) > component_tolerance
            or abs(components[1]) > component_tolerance
        ):
            raise ValueError(
                "NodalLoads must contain vertical Stage 02 self-weight only; "
                "index {0} has a horizontal component.".format(index)
            )
        weight = -components[2]
        if weight < -component_tolerance:
            raise ValueError(
                "NodalLoads contains a negative nodal weight at index {0}."
                .format(index)
            )
        weights.append(max(0.0, weight))
    if sum(weights) <= LOAD_COMPONENT_TOLERANCE_KN:
        raise ValueError("NodalLoads must contain positive total Stage 02 weight.")
    return weights


def _coerce_plane(value):
    if isinstance(value, rg.Plane):
        return rg.Plane(value)
    if hasattr(value, "Value") and isinstance(value.Value, rg.Plane):
        return rg.Plane(value.Value)
    return None


def _coerce_curve(value):
    if isinstance(value, rg.Curve):
        return value
    if hasattr(value, "Value") and isinstance(value.Value, rg.Curve):
        return value.Value
    return None


def _coerce_pair(value):
    if value is None:
        return None
    if hasattr(value, "Value"):
        return _coerce_pair(value.Value)
    try:
        items = list(value)
    except TypeError:
        return None
    if len(items) != 2:
        return None
    return str(items[0]), str(items[1])


def _finite_number(value, name, minimum=None):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError("{0} must be a number.".format(name))
    if not math.isfinite(number) or (
        minimum is not None and number < minimum
    ):
        raise ValueError(
            "{0} must be finite and at least {1}.".format(name, minimum)
        )
    return number


def _point_key(point, tolerance):
    scale = 1.0 / max(tolerance, 1.0e-9)
    return (
        int(round(point.X * scale)),
        int(round(point.Y * scale)),
        int(round(point.Z * scale)),
    )


def _ordered_network(lines, tolerance):
    points = []
    key_to_index = {}
    edges = []
    seen = set()

    def vertex_index(point):
        key = _point_key(point, tolerance)
        if key not in key_to_index:
            key_to_index[key] = len(points)
            points.append(rg.Point3d(point))
        return key_to_index[key]

    for line in lines:
        if not line.IsValid or line.Length <= tolerance:
            raise ValueError("TNAEdges contains an invalid or zero-length line.")
        a = vertex_index(line.From)
        b = vertex_index(line.To)
        key = tuple(sorted((a, b)))
        if key in seen:
            raise ValueError("TNAEdges contains a duplicate edge.")
        seen.add(key)
        edges.append((a, b))
    return points, edges


def _radial_topology(points, edges):
    adjacency = [[] for _ in points]
    for edge_index, (a, b) in enumerate(edges):
        adjacency[a].append((b, edge_index))
        adjacency[b].append((a, edge_index))
    supports = [
        index for index, neighbours in enumerate(adjacency)
        if len(neighbours) == 1 and points[index].Z <= 1.0e-6
    ]
    crowns = [
        index for index, neighbours in enumerate(adjacency)
        if len(neighbours) >= 3
    ]
    if len(crowns) != 1 or len(supports) < 3:
        raise ValueError("TNAEdges must form one radial crown-to-support tree.")
    if len(edges) != len(points) - 1:
        raise ValueError("TNAEdges must form a connected tree.")
    return adjacency, crowns[0], supports


def _recover_ideal_force_densities(points, edges, adjacency, supports, loads):
    support_set = set(supports)
    free = [index for index in range(len(points)) if index not in support_set]
    matrix = np.zeros((3 * len(free), len(edges)), dtype=float)
    target = np.zeros(3 * len(free), dtype=float)
    for local_index, node_index in enumerate(free):
        point = points[node_index]
        load = loads[node_index]
        target[3 * local_index:3 * local_index + 3] = (
            load.X, load.Y, load.Z
        )
        for neighbour, edge_index in adjacency[node_index]:
            other = points[neighbour]
            matrix[3 * local_index:3 * local_index + 3, edge_index] = (
                point.X - other.X,
                point.Y - other.Y,
                point.Z - other.Z,
            )
    q, residual_values, rank, singular_values = np.linalg.lstsq(
        matrix, target, rcond=None
    )
    residual = float(np.linalg.norm(matrix.dot(q) - target))
    total_weight = -sum(load.Z for load in loads)
    tolerance = max(1.0e-6, total_weight * 1.0e-6)
    forces = [
        -float(q[index]) * points[a].DistanceTo(points[b])
        for index, (a, b) in enumerate(edges)
    ]
    valid = (
        int(rank) == len(edges)
        and residual <= tolerance
        and all(float(value) < -COMPRESSION_EPSILON for value in q)
        and all(value > COMPRESSION_EPSILON for value in forces)
    )
    return q, forces, residual, tolerance, int(rank), free, valid


def _vector_equilibrium_matrix(points, edges, free):
    row_by_node = {node: index for index, node in enumerate(free)}
    matrix = np.zeros((3 * len(free), 3 * len(edges)), dtype=float)
    identity = np.eye(3)
    for edge_index, (a, b) in enumerate(edges):
        column = slice(3 * edge_index, 3 * edge_index + 3)
        if a in row_by_node:
            row = 3 * row_by_node[a]
            matrix[row:row + 3, column] += identity
        if b in row_by_node:
            row = 3 * row_by_node[b]
            matrix[row:row + 3, column] -= identity
    return matrix


def _ideal_edge_vectors(points, edges, force_densities):
    values = np.zeros(3 * len(edges), dtype=float)
    for edge_index, (a, b) in enumerate(edges):
        delta = points[a] - points[b]
        values[3 * edge_index:3 * edge_index + 3] = (
            float(force_densities[edge_index])
            * np.array([delta.X, delta.Y, delta.Z], dtype=float)
        )
    return values


def _profile_centroid(profile):
    properties = rg.AreaMassProperties.Compute(profile)
    if properties is None:
        raise RuntimeError("Could not compute a JointProfile centroid.")
    try:
        return rg.Point3d(properties.Centroid)
    finally:
        properties.Dispose()


def _line_plane_candidate(line, plane, tolerance):
    normal = rg.Vector3d(plane.Normal)
    direction = rg.Vector3d(line.Direction)
    if not normal.Unitize() or not direction.Unitize():
        return None
    denominator = rg.Vector3d.Multiply(line.Direction, normal)
    if abs(denominator) <= 1.0e-12:
        return None
    origin_delta = rg.Vector3d(
        plane.Origin.X - line.From.X,
        plane.Origin.Y - line.From.Y,
        plane.Origin.Z - line.From.Z,
    )
    parameter = rg.Vector3d.Multiply(origin_delta, normal) / denominator
    parameter_tolerance = tolerance / max(line.Length, tolerance)
    if parameter < -parameter_tolerance or parameter > 1.0 + parameter_tolerance:
        return None
    parameter = min(max(parameter, 0.0), 1.0)
    point = line.PointAt(parameter)
    alignment = abs(rg.Vector3d.Multiply(direction, normal))
    alignment = min(max(alignment, -1.0), 1.0)
    return point, math.degrees(math.acos(alignment))


def _map_joint(profile, plane, lines, tolerance):
    centroid = _profile_centroid(profile)
    candidates = []
    for edge_index, line in enumerate(lines):
        candidate = _line_plane_candidate(line, plane, tolerance)
        if candidate is None:
            continue
        point, angle_error = candidate
        if angle_error > MAX_ALIGNMENT_ERROR_DEGREES:
            continue
        containment = profile.Contains(point, plane, tolerance)
        if containment not in (
            rg.PointContainment.Inside,
            rg.PointContainment.Coincident,
        ):
            continue
        candidates.append(
            (angle_error, point.DistanceTo(centroid), edge_index, point)
        )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    return candidates[0]


def _free_load_target(loads, free):
    target = np.zeros(3 * len(free), dtype=float)
    for local_index, node_index in enumerate(free):
        load = loads[node_index]
        target[3 * local_index:3 * local_index + 3] = (
            load.X, load.Y, load.Z
        )
    return target


def _horizontal_load_targets(loads, free):
    target_x = np.zeros(3 * len(free), dtype=float)
    target_y = np.zeros(3 * len(free), dtype=float)
    for local_index, node_index in enumerate(free):
        weight = max(0.0, -loads[node_index].Z)
        target_x[3 * local_index] = weight
        target_y[3 * local_index + 1] = weight
    return target_x, target_y


def _prepare_joint_resultants(
    tna_edges,
    nodal_loads,
    joint_planes,
    joint_profiles,
    joint_ids,
    joint_block_pairs,
):
    tolerance = _document_tolerance()
    lines = [_coerce_line(value) for value in _as_list(tna_edges)]
    loads = [_coerce_vector(value) for value in _as_list(nodal_loads)]
    planes = [_coerce_plane(value) for value in _as_list(joint_planes)]
    profiles = [_coerce_curve(value) for value in _as_list(joint_profiles)]
    identifiers = [str(value) for value in _as_list(joint_ids)]
    pairs = [_coerce_pair(value) for value in _as_list(joint_block_pairs)]
    if not lines or any(value is None for value in lines):
        raise ValueError("TNAEdges must contain valid lines.")
    if any(value is None for value in loads):
        raise ValueError("NodalLoads must contain valid vectors.")
    joint_count = len(planes)
    if (
        not joint_count
        or len(profiles) != joint_count
        or len(identifiers) != joint_count
        or len(pairs) != joint_count
        or any(value is None for value in planes + profiles)
        or any(value is None or value[0] == value[1] for value in pairs)
    ):
        raise ValueError(
            "JointPlanes, JointProfiles, JointIDs and JointBlockPairs must "
            "correspond one-to-one."
        )
    if len(set(identifiers)) != joint_count:
        raise ValueError("JointIDs must be unique.")

    points, edges = _ordered_network(
        lines, max(tolerance * 0.1, 1.0e-7)
    )
    if len(loads) != len(points):
        raise ValueError("NodalLoads must correspond to the ordered TNA nodes.")
    nodal_weights = _validated_stage02_weights(loads)
    adjacency, crown, supports = _radial_topology(points, edges)
    (
        force_densities,
        ideal_forces,
        ideal_residual,
        ideal_tolerance,
        ideal_rank,
        free,
        ideal_valid,
    ) = _recover_ideal_force_densities(
        points, edges, adjacency, supports, loads
    )
    total_nodal_weight = sum(nodal_weights)
    free_node_weight = sum(nodal_weights[index] for index in free)
    support_node_weight = sum(nodal_weights[index] for index in supports)
    weight_partition_error = abs(
        total_nodal_weight - free_node_weight - support_node_weight
    )
    weight_partition_tolerance = max(
        LOAD_COMPONENT_TOLERANCE_KN,
        total_nodal_weight * WEIGHT_PARTITION_RELATIVE_TOLERANCE,
    )
    if weight_partition_error > weight_partition_tolerance:
        raise RuntimeError(
            "Stage 02 nodal weights do not close between free and support nodes."
        )

    contact_points = []
    angle_errors = []
    mapped_edges = []
    failures = []
    for index, (plane, profile) in enumerate(zip(planes, profiles)):
        mapped = _map_joint(profile, plane, lines, tolerance)
        if mapped is None:
            failures.append((
                index, identifiers[index], "no_TNA_joint_intersection"
            ))
            continue
        angle_error, centroid_distance, edge_index, point = mapped
        contact_points.append(point)
        angle_errors.append(angle_error)
        mapped_edges.append(edge_index)

    mapping_valid = not failures and len(contact_points) == joint_count
    equilibrium = None
    ideal_vectors = None
    correction_basis = None
    vector_rank = 0
    if mapping_valid:
        equilibrium = _vector_equilibrium_matrix(points, edges, free)
        ideal_vectors = _ideal_edge_vectors(points, edges, force_densities)
        base_target = _free_load_target(loads, free)
        horizontal_x, horizontal_y = _horizontal_load_targets(loads, free)
        correction_targets = np.column_stack((
            base_target - equilibrium.dot(ideal_vectors),
            horizontal_x,
            horizontal_y,
        ))
        correction_basis, residual_values, vector_rank, singular_values = (
            np.linalg.lstsq(equilibrium, correction_targets, rcond=None)
        )

    return {
        "tolerance": tolerance,
        "lines": lines,
        "loads": loads,
        "planes": planes,
        "profiles": profiles,
        "identifiers": identifiers,
        "pairs": pairs,
        "joint_count": joint_count,
        "points": points,
        "edges": edges,
        "crown": crown,
        "supports": supports,
        "free": free,
        "force_densities": force_densities,
        "ideal_forces": ideal_forces,
        "ideal_residual": ideal_residual,
        "ideal_tolerance": ideal_tolerance,
        "ideal_rank": ideal_rank,
        "ideal_valid": ideal_valid,
        "contact_points": contact_points,
        "angle_errors": angle_errors,
        "mapped_edges": mapped_edges,
        "mapping_failures": failures,
        "mapping_valid": mapping_valid,
        "equilibrium": equilibrium,
        "ideal_vectors": ideal_vectors,
        "correction_basis": correction_basis,
        "vector_rank": int(vector_rank),
        "total_nodal_weight": total_nodal_weight,
        "free_node_weight": free_node_weight,
        "support_node_weight": support_node_weight,
        "weight_partition_error": weight_partition_error,
        "weight_partition_tolerance": weight_partition_tolerance,
    }


def _evaluate_prepared_joint_resultants(
    prepared,
    horizontal_load_ratio,
    horizontal_load_angle,
):
    ratio = _finite_number(
        horizontal_load_ratio, "HorizontalLoadRatio", minimum=0.0
    )
    angle = _finite_number(
        horizontal_load_angle, "HorizontalLoadAngle"
    )
    radians = math.radians(angle)
    direction = rg.Vector3d(math.cos(radians), math.sin(radians), 0.0)

    lines = prepared["lines"]
    loads = prepared["loads"]
    planes = prepared["planes"]
    identifiers = prepared["identifiers"]
    joint_count = prepared["joint_count"]
    points = prepared["points"]
    edges = prepared["edges"]
    crown = prepared["crown"]
    supports = prepared["supports"]
    free = prepared["free"]
    ideal_forces = prepared["ideal_forces"]
    ideal_residual = prepared["ideal_residual"]
    ideal_tolerance = prepared["ideal_tolerance"]
    ideal_rank = prepared["ideal_rank"]
    ideal_valid = prepared["ideal_valid"]
    contact_points = prepared["contact_points"]
    angle_errors = prepared["angle_errors"]
    mapped_edges = prepared["mapped_edges"]
    failures = list(prepared["mapping_failures"])
    mapping_valid = prepared["mapping_valid"]

    contact_resultants = []
    tangential_forces = []
    normal_forces = []
    load_case = "ideal" if ratio == 0.0 else "horizontal_sensitivity"
    equilibrium_valid = True
    vector_rank_text = "not_applicable"
    residual = ideal_residual
    residual_tolerance = ideal_tolerance
    correction_norm = 0.0
    total_nodal_weight = prepared["total_nodal_weight"]
    free_node_weight = prepared["free_node_weight"]
    support_node_weight = prepared["support_node_weight"]
    free_horizontal_resultant = ratio * free_node_weight
    full_structure_inertia_resultant = ratio * total_nodal_weight
    support_local_inertia_resultant = ratio * support_node_weight
    ground_acceleration_angle = (angle + 180.0) % 360.0

    if mapping_valid and ratio == 0.0:
        for plane, edge_index in zip(planes, mapped_edges):
            normal = rg.Vector3d(plane.Normal)
            normal.Unitize()
            magnitude = ideal_forces[edge_index]
            contact_resultants.append(normal * magnitude)
            tangential_forces.append(0.0)
            normal_forces.append(magnitude)
    elif mapping_valid:
        equilibrium = prepared["equilibrium"]
        ideal_vectors = prepared["ideal_vectors"]
        correction_basis = prepared["correction_basis"]
        correction = (
            correction_basis[:, 0]
            + ratio * (
                math.cos(radians) * correction_basis[:, 1]
                + math.sin(radians) * correction_basis[:, 2]
            )
        )
        solved_vectors = ideal_vectors + correction
        target = _free_load_target(loads, free)
        horizontal_x, horizontal_y = _horizontal_load_targets(loads, free)
        target += ratio * (
            math.cos(radians) * horizontal_x
            + math.sin(radians) * horizontal_y
        )
        residual = float(np.linalg.norm(
            equilibrium.dot(solved_vectors) - target
        ))
        total_weight = -sum(load.Z for load in loads)
        residual_tolerance = max(1.0e-7, total_weight * 1.0e-8)
        equilibrium_valid = (
            prepared["vector_rank"] == equilibrium.shape[0]
            and residual <= residual_tolerance
        )
        vector_rank_text = "{0}/{1}".format(
            prepared["vector_rank"], equilibrium.shape[0]
        )
        correction_norm = float(np.linalg.norm(correction))
        for index, (plane, edge_index) in enumerate(zip(planes, mapped_edges)):
            vector = rg.Vector3d(
                solved_vectors[3 * edge_index],
                solved_vectors[3 * edge_index + 1],
                solved_vectors[3 * edge_index + 2],
            )
            normal = rg.Vector3d(plane.Normal)
            normal.Unitize()
            if rg.Vector3d.Multiply(vector, normal) < 0.0:
                vector = -vector
            normal_force = rg.Vector3d.Multiply(vector, normal)
            tangential = vector - normal * normal_force
            if normal_force <= COMPRESSION_EPSILON:
                failures.append((
                    index, identifiers[index], "nonpositive_normal_force"
                ))
                continue
            contact_resultants.append(vector)
            tangential_forces.append(tangential.Length)
            normal_forces.append(normal_force)

    mapping_valid = (
        mapping_valid
        and not failures
        and len(contact_resultants) == joint_count
        and len(tangential_forces) == joint_count
    )
    valid = ideal_valid and equilibrium_valid and mapping_valid
    report = (
        f"stage=joint_resultants; load_case={load_case}; nodes={len(points)}; "
        f"edges={len(edges)}; free_nodes={len(free)}; supports={len(supports)}; "
        f"crown_node={crown}; horizontal_load_ratio={ratio:.6g}; "
        "kh_definition=a_h_over_g; "
        f"equivalent_inertia_force_angle_deg={angle:.6g}; "
        f"corresponding_base_acceleration_angle_deg="
        f"{ground_acceleration_angle:.6g}; "
        f"horizontal_direction=({direction.X:.6g},{direction.Y:.6g},0); "
        f"stage02_total_nodal_weight_kN={total_nodal_weight:.12g}; "
        f"free_node_weight_kN={free_node_weight:.12g}; "
        f"support_node_weight_kN={support_node_weight:.12g}; "
        f"weight_partition_error_kN="
        f"{prepared['weight_partition_error']:.12g}; "
        f"free_node_horizontal_target_resultant_kN="
        f"{free_horizontal_resultant:.12g}; "
        f"full_structure_equivalent_inertia_resultant_kN="
        f"{full_structure_inertia_resultant:.12g}; "
        f"support_local_inertia_resultant_kN="
        f"{support_local_inertia_resultant:.12g}; "
        f"ideal_force_density_rank={ideal_rank}/{len(edges)}; "
        f"ideal_equilibrium_residual_kN={ideal_residual:.6g}; "
        f"ideal_valid={ideal_valid}; "
        f"vector_equilibrium_rank={vector_rank_text}; "
        f"vector_equilibrium_residual_kN={residual:.6g}; "
        f"vector_equilibrium_tolerance_kN={residual_tolerance:.6g}; "
        f"force_vector_correction_norm_kN={correction_norm:.6g}; "
        f"joints_mapped={len(contact_resultants)}/{joint_count}; "
        f"normal_force_min_max_kN="
        f"({min(normal_forces) if normal_forces else 0.0:.6g},"
        f"{max(normal_forces) if normal_forces else 0.0:.6g}); "
        f"tangential_force_min_max_kN="
        f"({min(tangential_forces) if tangential_forces else 0.0:.6g},"
        f"{max(tangential_forces) if tangential_forces else 0.0:.6g}); "
        f"max_alignment_error_deg="
        f"{max(angle_errors) if angle_errors else 0.0:.6g}; "
        f"mapped_tna_edge_indices={mapped_edges}; failures={failures}; "
        f"equilibrium_valid={equilibrium_valid}; "
        f"mapping_valid={mapping_valid}; valid={valid}; "
        f"solution={'exact_ideal_axial' if ratio == 0.0 else 'min_change_from_ideal_axial_vectors'}; "
        "resultant_convention=second_block_on_first_block; "
        "case_angle_convention=equivalent_inertia_force_direction; "
        "base_acceleration_direction=(inertia_force_angle+180deg)_mod_360; "
        "support_node_inertia_treatment=local_support_reaction_not_free_node_solve; "
        f"scope={'TNA_compression_mapping_only' if ratio == 0.0 else 'mass_proportional_uniform_base_acceleration_full_direction_equivalent_static_inertia_sensitivity'}; "
        "free_TNA_node_translational_equilibrium_only; "
        "block_moments_support_contact_distribution_and_dynamic_response_not_checked; "
        "not_earthquake_safety_response_spectrum_or_time_history_analysis"
    )
    return (
        list(contact_points),
        contact_resultants,
        tangential_forces,
        valid,
        report,
    )


def create_joint_resultants(
    tna_edges,
    nodal_loads,
    joint_planes,
    joint_profiles,
    joint_ids,
    joint_block_pairs,
    horizontal_load_ratio,
    horizontal_load_angle,
):
    prepared = _prepare_joint_resultants(
        tna_edges,
        nodal_loads,
        joint_planes,
        joint_profiles,
        joint_ids,
        joint_block_pairs,
    )
    return _evaluate_prepared_joint_resultants(
        prepared,
        horizontal_load_ratio,
        horizontal_load_angle,
    )


def _direction_angles(horizontal_load_ratio, angle_step_degrees):
    ratio = _finite_number(
        horizontal_load_ratio, "HorizontalLoadRatio", minimum=0.0
    )
    if ratio == 0.0:
        return ratio, [0.0]
    step = _finite_number(
        angle_step_degrees, "HorizontalAngleStepDegrees", minimum=1.0e-9
    )
    divisions = 360.0 / step
    rounded = int(round(divisions))
    if rounded < 2 or rounded > 720 or abs(divisions - rounded) > 1.0e-9:
        raise ValueError(
            "HorizontalAngleStepDegrees must divide 360 exactly and produce "
            "between 2 and 720 directional cases."
        )
    return ratio, [index * step for index in range(rounded)]


def _canonical_angle(angle):
    value = float(angle) % 360.0
    if abs(value - 360.0) <= 1.0e-9 or abs(value) <= 1.0e-9:
        value = 0.0
    return round(value, 10)


def _refinement_parameters(
    coarse_step_degrees,
    refinement_step_degrees,
    refinement_half_width_degrees,
    friction_coefficient,
    sliding_safety_factor,
):
    coarse_step = _finite_number(
        coarse_step_degrees, "HorizontalAngleStepDegrees", minimum=1.0e-9
    )
    refinement_step = _finite_number(
        refinement_step_degrees,
        "RefinementAngleStepDegrees",
        minimum=1.0e-9,
    )
    half_width = _finite_number(
        refinement_half_width_degrees,
        "RefinementHalfWidthDegrees",
        minimum=refinement_step,
    )
    friction = _finite_number(
        friction_coefficient, "FrictionCoefficient", minimum=1.0e-12
    )
    sliding_factor = _finite_number(
        sliding_safety_factor, "SlidingSafetyFactor", minimum=1.0e-12
    )
    coarse_ratio = coarse_step / refinement_step
    half_width_ratio = half_width / refinement_step
    full_circle_ratio = 360.0 / refinement_step
    if (
        abs(coarse_ratio - round(coarse_ratio)) > 1.0e-9
        or abs(half_width_ratio - round(half_width_ratio)) > 1.0e-9
        or abs(full_circle_ratio - round(full_circle_ratio)) > 1.0e-9
    ):
        raise ValueError(
            "RefinementAngleStepDegrees must divide the coarse step, "
            "RefinementHalfWidthDegrees and 360 exactly."
        )
    if refinement_step > coarse_step + 1.0e-12:
        raise ValueError(
            "RefinementAngleStepDegrees must not exceed "
            "HorizontalAngleStepDegrees."
        )
    if int(round(full_circle_ratio)) > 720:
        raise ValueError(
            "RefinementAngleStepDegrees must produce no more than 720 "
            "directions over 360 degrees."
        )
    return refinement_step, half_width, friction, sliding_factor


def _signed_sliding_deficit(resultant, plane, friction, sliding_factor):
    normal = rg.Vector3d(plane.Normal)
    if not normal.Unitize():
        raise ValueError("JointPlanes contains a degenerate plane.")
    force = rg.Vector3d(resultant)
    normal_force = rg.Vector3d.Multiply(force, normal)
    tangent = force - normal * normal_force
    return sliding_factor * tangent.Length - friction * normal_force


def _adaptive_refinement_window(
    coarse_peak_angle,
    score_at_angle,
    refinement_step_degrees,
    refinement_half_width_degrees,
):
    step = float(refinement_step_degrees)
    half_width = float(refinement_half_width_degrees)
    lower = float(coarse_peak_angle) - half_width
    upper = float(coarse_peak_angle) + half_width
    evaluated = set()
    expansions = 0
    full_circle_covered = False

    while True:
        count = int(round((upper - lower) / step))
        unwrapped = [lower + index * step for index in range(count + 1)]
        for angle in unwrapped:
            evaluated.add(_canonical_angle(angle))
        scores = [float(score_at_angle(angle)) for angle in unwrapped]
        maximum = max(scores)
        lower_is_peak = (
            scores[0] >= maximum - REFINEMENT_SCORE_TOLERANCE_KN
        )
        upper_is_peak = (
            scores[-1] >= maximum - REFINEMENT_SCORE_TOLERANCE_KN
        )
        if not lower_is_peak and not upper_is_peak:
            break
        if upper - lower + step >= 360.0 - 1.0e-9:
            full_circle_covered = True
            break
        if lower_is_peak:
            lower -= half_width
        if upper_is_peak:
            upper += half_width
        expansions += 1

    final_count = int(round((upper - lower) / step))
    final_unwrapped = [lower + index * step for index in range(final_count + 1)]
    final_scores = [float(score_at_angle(angle)) for angle in final_unwrapped]
    maximum = max(final_scores)
    peak_candidates = [
        angle
        for angle, score in zip(final_unwrapped, final_scores)
        if score >= maximum - REFINEMENT_SCORE_TOLERANCE_KN
    ]
    interior_candidates = [
        angle
        for angle in peak_candidates
        if angle > lower + 1.0e-9 and angle < upper - 1.0e-9
    ]
    selected = (
        min(interior_candidates, key=lambda value: _canonical_angle(value))
        if interior_candidates
        else min(peak_candidates, key=lambda value: _canonical_angle(value))
    )
    return (
        sorted(evaluated),
        _canonical_angle(selected),
        maximum,
        expansions,
        full_circle_covered,
    )


def create_directional_joint_resultants(
    tna_edges,
    nodal_loads,
    joint_planes,
    joint_profiles,
    joint_ids,
    joint_block_pairs,
    horizontal_load_ratio,
    horizontal_angle_step_degrees,
    refinement_angle_step_degrees=1.0,
    refinement_half_width_degrees=5.0,
    friction_coefficient=0.55,
    sliding_safety_factor=1.5,
):
    ratio, coarse_angles = _direction_angles(
        horizontal_load_ratio, horizontal_angle_step_degrees
    )
    identifiers = [str(value) for value in _as_list(joint_ids)]
    planes = [_coerce_plane(value) for value in _as_list(joint_planes)]
    if (
        not identifiers
        or len(planes) != len(identifiers)
        or any(value is None for value in planes)
    ):
        raise ValueError("JointPlanes and JointIDs must correspond one-to-one.")
    if ratio > 0.0:
        (
            refinement_step,
            refinement_half_width,
            friction,
            sliding_factor,
        ) = _refinement_parameters(
            horizontal_angle_step_degrees,
            refinement_angle_step_degrees,
            refinement_half_width_degrees,
            friction_coefficient,
            sliding_safety_factor,
        )
    else:
        refinement_step = 0.0
        refinement_half_width = 0.0
        friction = 0.0
        sliding_factor = 0.0

    result_cache = {}
    case_failures = []
    prepared = _prepare_joint_resultants(
        tna_edges,
        nodal_loads,
        joint_planes,
        joint_profiles,
        joint_ids,
        joint_block_pairs,
    )

    def evaluate(angle):
        canonical = _canonical_angle(angle)
        if canonical not in result_cache:
            result = _evaluate_prepared_joint_resultants(
                prepared,
                ratio,
                canonical,
            )
            result_cache[canonical] = result
            if not result[3]:
                case_failures.append((canonical, result[4]))
        return result_cache[canonical]

    for angle in coarse_angles:
        evaluate(angle)

    coarse_peak_angles = []
    coarse_candidate_angles = []
    refined_peak_angles = []
    refined_peak_deficits = []
    boundary_expansions = []
    full_circle_refinements = []
    if ratio > 0.0 and not case_failures:
        for joint_index, plane in enumerate(planes):
            coarse_scores = [
                _signed_sliding_deficit(
                    evaluate(angle)[1][joint_index],
                    plane,
                    friction,
                    sliding_factor,
                )
                for angle in coarse_angles
            ]
            coarse_maximum = max(coarse_scores)
            coarse_peak = min(
                angle
                for angle, score in zip(coarse_angles, coarse_scores)
                if score >= coarse_maximum - REFINEMENT_SCORE_TOLERANCE_KN
            )
            coarse_peak_angles.append(coarse_peak)
            candidate_indices = []
            for index, score in enumerate(coarse_scores):
                previous = coarse_scores[index - 1]
                following = coarse_scores[(index + 1) % len(coarse_scores)]
                is_not_lower = (
                    score >= previous - REFINEMENT_SCORE_TOLERANCE_KN
                    and score >= following - REFINEMENT_SCORE_TOLERANCE_KN
                )
                is_not_flat = (
                    score > previous + REFINEMENT_SCORE_TOLERANCE_KN
                    or score > following + REFINEMENT_SCORE_TOLERANCE_KN
                )
                if is_not_lower and is_not_flat:
                    candidate_indices.append(index)
            if not candidate_indices:
                candidate_indices = [coarse_angles.index(coarse_peak)]
            candidate_angles = [
                coarse_angles[index] for index in candidate_indices
            ]
            if coarse_peak not in candidate_angles:
                candidate_angles.append(coarse_peak)
                candidate_angles.sort()
            coarse_candidate_angles.append(candidate_angles)

            def score_at_angle(angle, index=joint_index, joint_plane=plane):
                result = evaluate(angle)
                if not result[3]:
                    return float("-inf")
                return _signed_sliding_deficit(
                    result[1][index],
                    joint_plane,
                    friction,
                    sliding_factor,
                )

            candidate_results = []
            for candidate_angle in candidate_angles:
                result = _adaptive_refinement_window(
                    candidate_angle,
                    score_at_angle,
                    refinement_step,
                    refinement_half_width,
                )
                for angle in result[0]:
                    evaluate(angle)
                candidate_results.append(result)
            best_deficit = max(value[2] for value in candidate_results)
            best_candidates = [
                value
                for value in candidate_results
                if value[2]
                >= best_deficit - REFINEMENT_SCORE_TOLERANCE_KN
            ]
            best = min(best_candidates, key=lambda value: value[1])
            refined_peak_angles.append(best[1])
            refined_peak_deficits.append(best[2])
            boundary_expansions.append(
                sum(value[3] for value in candidate_results)
            )
            full_circle_refinements.append(
                any(value[4] for value in candidate_results)
            )
    elif ratio == 0.0:
        coarse_peak_angles = [0.0] * len(identifiers)
        coarse_candidate_angles = [[0.0] for value in identifiers]
        refined_peak_angles = [0.0] * len(identifiers)
        refined_peak_deficits = [0.0] * len(identifiers)
        boundary_expansions = [0] * len(identifiers)
        full_circle_refinements = [False] * len(identifiers)

    angles = sorted(result_cache)
    case_joint_ids = []
    contact_points = []
    contact_resultants = []
    tangential_forces = []
    for angle in angles:
        result = result_cache[angle]
        if not result[3]:
            continue
        case_joint_ids.extend(identifiers)
        contact_points.extend(result[0])
        contact_resultants.extend(result[1])
        tangential_forces.extend(result[2])

    expected_rows = len(angles) * len(identifiers)
    total_nodal_weight = prepared["total_nodal_weight"]
    free_node_weight = prepared["free_node_weight"]
    support_node_weight = prepared["support_node_weight"]
    valid = (
        not case_failures
        and len(case_joint_ids)
        == len(contact_points)
        == len(contact_resultants)
        == len(tangential_forces)
        == expected_rows
    )
    report = (
        f"stage=directional_joint_resultants; horizontal_load_ratio={ratio:.6g}; "
        "kh_definition=a_h_over_g; "
        "case_angle_convention=equivalent_horizontal_inertia_force_direction; "
        "base_acceleration_direction=(CaseAngle+180deg)_mod_360; "
        f"stage02_total_nodal_weight_kN={total_nodal_weight:.12g}; "
        f"free_node_weight_kN={free_node_weight:.12g}; "
        f"support_node_weight_kN={support_node_weight:.12g}; "
        f"weight_partition_error_kN="
        f"{prepared['weight_partition_error']:.12g}; "
        f"free_node_horizontal_target_resultant_kN="
        f"{ratio * free_node_weight:.12g}; "
        f"full_structure_equivalent_inertia_resultant_kN="
        f"{ratio * total_nodal_weight:.12g}; "
        f"support_local_inertia_resultant_kN="
        f"{ratio * support_node_weight:.12g}; "
        "support_node_inertia_treatment=local_support_reaction_not_free_node_solve; "
        f"directional_sweep={ratio > 0.0}; angle_cases={len(angles)}; "
        f"coarse_angle_cases={len(coarse_angles)}; coarse_angle_step_deg="
        f"{horizontal_angle_step_degrees if ratio > 0.0 else 0.0}; "
        f"refinement_step_deg={refinement_step:.6g}; "
        f"refinement_half_width_deg={refinement_half_width:.6g}; "
        f"angle_min_max_deg=({min(angles):.6g},{max(angles):.6g}); "
        f"joints_per_case={len(identifiers)}; flattened_rows={len(case_joint_ids)}; "
        f"coarse_peak_angles_deg={coarse_peak_angles}; "
        f"coarse_local_peak_candidates_deg={coarse_candidate_angles}; "
        f"refined_peak_angles_deg={refined_peak_angles}; "
        f"refined_peak_signed_deficits_kN={refined_peak_deficits}; "
        f"boundary_expansions={boundary_expansions}; "
        f"full_circle_refinements={full_circle_refinements}; "
        f"case_failures={case_failures}; valid={valid}; "
        "adaptive_refinement=all_coarse_local_maxima_of_per_joint_final_sliding_deficit; "
        "boundary_rule=extend_while_boundary_equals_current_maximum; "
        "flattening=case_major_each_case_follows_jointids_order; "
        "equilibrium_scope=free_TNA_node_translational_force_balance_only; "
        "block_moments_support_contact_distribution_and_dynamic_response_not_checked; "
        "scope=mass_proportional_uniform_base_acceleration_full_direction_equivalent_static_inertia_sensitivity_not_earthquake_safety_check"
    )
    return (
        angles,
        case_joint_ids,
        contact_points,
        contact_resultants,
        tangential_forces,
        valid,
        report,
    )


try:
    _component_inputs = (
        TNAEdges,
        NodalLoads,
        JointPlanes,
        JointProfiles,
        JointIDs,
        JointBlockPairs,
        HorizontalLoadRatio,
        HorizontalAngleStepDegrees,
        RefinementAngleStepDegrees,
        RefinementHalfWidthDegrees,
        FrictionCoefficient,
        SlidingSafetyFactor,
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        CaseAngles,
        CaseJointIDs,
        ContactPoints,
        ContactResultants,
        TangentialForces,
        Valid,
        Report,
    ) = create_directional_joint_resultants(*_component_inputs)
