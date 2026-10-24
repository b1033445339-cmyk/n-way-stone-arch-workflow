#! python 3
"""Stage 07F - per-joint sliding-deficit envelope (paper version).

Inputs: JointPlanes, JointIDs, CaseAngles, CaseJointIDs, ContactResultants,
PressureCentroids, SlidingCaseKeys, FrictionCoefficient, SlidingSafetyFactor.
Outputs: DeficitJointIDs, DeficitVectors, SlidingDeficits,
GoverningHorizontalAngles, GoverningNormalForces, GoverningTangentialForces,
GoverningSlidingUtilisations, GoverningPressureCentroids, Valid, Report.

D = max(0, SlidingSafetyFactor*T - FrictionCoefficient*N), in kN.
Each joint retains D and the paired quantities from its own governing case.
Different joints need not share a governing direction. Positive D identifies
an interface for review, not connector demand/capacity or structural risk.
The vector opposes the tangential resultant on the first adjacent block.
The numerical kernel and tie-breaking rules are retained from the verified
07F implementation. No connector sizing or contact-strength verdict is made.
"""

import math

import Rhino.Geometry as rg


DEMAND_TOLERANCE_KN = 1.0e-9


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]


def _coerce_plane(value):
    if isinstance(value, rg.Plane):
        return rg.Plane(value)
    if hasattr(value, "Value") and isinstance(value.Value, rg.Plane):
        return rg.Plane(value.Value)
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


def _coerce_point(value):
    if isinstance(value, rg.Point3d):
        return rg.Point3d(value)
    if hasattr(value, "Value") and isinstance(value.Value, rg.Point3d):
        return rg.Point3d(value.Value)
    return None


def _positive_number(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError("{0} must be a number.".format(name))
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("{0} must be finite and greater than zero.".format(name))
    return number


def create_sliding_deficit(
    joint_planes,
    joint_ids,
    contact_resultants,
    sliding_joint_ids,
    friction_coefficient,
    sliding_safety_factor,
    _prepared_normals=None,
):
    planes = [_coerce_plane(value) for value in _as_list(joint_planes)]
    identifiers = [str(value) for value in _as_list(joint_ids)]
    resultants = [_coerce_vector(value) for value in _as_list(contact_resultants)]
    supplied_sliding_ids = [
        str(value) for value in _as_list(sliding_joint_ids)
    ]

    count = len(planes)
    if not count or not all(
        len(values) == count for values in (identifiers, resultants)
    ):
        raise ValueError("All Stage 07F joint inputs must correspond one-to-one.")
    if any(value is None for value in planes + resultants):
        raise ValueError("Stage 07F contains invalid Plane or Vector data.")
    if len(set(identifiers)) != count:
        raise ValueError("JointIDs must be unique.")
    if len(set(supplied_sliding_ids)) != len(supplied_sliding_ids):
        raise ValueError("SlidingJointIDs must be unique.")
    unknown_ids = [
        value for value in supplied_sliding_ids if value not in set(identifiers)
    ]
    if unknown_ids:
        raise ValueError(
            "SlidingJointIDs contains unknown IDs: {0}".format(unknown_ids)
        )
    if _prepared_normals is None:
        prepared_normals = [None] * count
    else:
        prepared_normals = list(_prepared_normals)
        if len(prepared_normals) != count:
            raise ValueError("Prepared Stage 07F normals must match JointIDs.")

    friction = _positive_number(friction_coefficient, "FrictionCoefficient")
    sliding_factor = _positive_number(
        sliding_safety_factor, "SlidingSafetyFactor"
    )

    deficit_ids = []
    demand_vectors = []
    demand_forces = []
    recomputed_sliding_ids = []
    normal_forces = []
    tangential_forces = []
    friction_resistances = []
    failures = []

    for index, (plane, resultant, prepared_normal) in enumerate(zip(
        planes, resultants, prepared_normals
    )):
        try:
            if prepared_normal is None:
                normal = rg.Vector3d(plane.Normal)
                if not normal.Unitize():
                    raise ValueError("JointPlane has an invalid normal.")
            else:
                normal = prepared_normal
            normal_force = rg.Vector3d.Multiply(resultant, normal)
            if not math.isfinite(normal_force) or normal_force <= 0.0:
                raise ValueError(
                    "ContactResultant has no positive compressive normal force."
                )
            tangential = resultant - normal * normal_force
            tangential_force = tangential.Length
            if not math.isfinite(tangential_force):
                raise ValueError("ContactResultant has an invalid tangent force.")

            friction_resistance = friction * normal_force
            demand = max(
                0.0,
                sliding_factor * tangential_force - friction_resistance,
            )
            normal_forces.append(normal_force)
            tangential_forces.append(tangential_force)
            friction_resistances.append(friction_resistance)
            if demand > DEMAND_TOLERANCE_KN:
                direction = rg.Vector3d(tangential)
                if not direction.Unitize():
                    raise ValueError(
                        "A nonzero sliding deficit has no valid resisting direction."
                    )
                identifier = identifiers[index]
                recomputed_sliding_ids.append(identifier)
                deficit_ids.append(identifier)
                demand_vectors.append(-direction * demand)
                demand_forces.append(demand)
        except Exception as error:
            failures.append((index, identifiers[index], str(error)))

    mapping_valid = recomputed_sliding_ids == supplied_sliding_ids
    if not mapping_valid:
        failures.append((
            "sliding_id_mismatch",
            "supplied={0}; recomputed={1}".format(
                supplied_sliding_ids, recomputed_sliding_ids
            ),
        ))
    valid = (
        not failures
        and mapping_valid
        and len(deficit_ids) == len(demand_vectors) == len(demand_forces)
    )
    report = (
        f"stage=deficit_demand; joints={count}; "
        "parameters_source=assumed_not_measured; "
        f"friction_coefficient={friction:.6g}; "
        f"sliding_safety_factor={sliding_factor:.6g}; "
        f"supplied_sliding_joints={len(supplied_sliding_ids)}; "
        f"deficit_joints={len(deficit_ids)}; deficit_joint_ids={deficit_ids}; "
        f"deficit_demand_force_min_max_kN="
        f"({min(demand_forces) if demand_forces else 0.0:.6g},"
        f"{max(demand_forces) if demand_forces else 0.0:.6g}); "
        f"normal_force_min_max_kN="
        f"({min(normal_forces) if normal_forces else 0.0:.6g},"
        f"{max(normal_forces) if normal_forces else 0.0:.6g}); "
        f"max_tangential_force_kN="
        f"{max(tangential_forces) if tangential_forces else 0.0:.6g}; "
        f"friction_resistance_min_max_kN="
        f"({min(friction_resistances) if friction_resistances else 0.0:.6g},"
        f"{max(friction_resistances) if friction_resistances else 0.0:.6g}); "
        f"mapping_valid={mapping_valid}; failures={failures}; valid={valid}; "
        "demand_model=max(0,safety_factor*tangent-friction*normal); "
        "deficit_vector_convention=resistance_on_first_block_opposes_tangent; "
        "pressure_is_assessed_separately; "
        "no_connector_capacity_assessed; input_load_case_not_independently_validated"
    )
    return deficit_ids, demand_vectors, demand_forces, valid, report


def _case_key(case_index, angle, identifier):
    return "case[{0}]/angle[{1:.6g}]/{2}".format(
        case_index, angle, identifier
    )


def create_directional_sliding_deficit_envelope(
    joint_planes,
    joint_ids,
    case_angles,
    case_joint_ids,
    contact_resultants,
    pressure_centroids,
    sliding_case_keys,
    friction_coefficient,
    sliding_safety_factor,
):
    planes = [_coerce_plane(value) for value in _as_list(joint_planes)]
    identifiers = [str(value) for value in _as_list(joint_ids)]
    try:
        angles = [float(value) for value in _as_list(case_angles)]
    except (TypeError, ValueError):
        raise ValueError("CaseAngles must contain numeric degree values.")
    row_ids = [str(value) for value in _as_list(case_joint_ids)]
    resultants = [_coerce_vector(value) for value in _as_list(contact_resultants)]
    centroids = [_coerce_point(value) for value in _as_list(pressure_centroids)]
    supplied_keys = [str(value) for value in _as_list(sliding_case_keys)]
    joint_count = len(identifiers)
    expected_rows = len(angles) * joint_count
    if (
        not joint_count
        or not angles
        or len(planes) != joint_count
        or len(set(identifiers)) != joint_count
        or len(row_ids) != expected_rows
        or len(resultants) != expected_rows
        or len(centroids) != expected_rows
        or any(value is None for value in planes + resultants + centroids)
    ):
        raise ValueError(
            "Stage 07F directional inputs must contain one complete "
            "JointIDs-sized block for every CaseAngles item."
        )
    if len(set(supplied_keys)) != len(supplied_keys):
        raise ValueError("SlidingCaseKeys must be unique.")
    identifier_indices = {
        identifier: index for index, identifier in enumerate(identifiers)
    }
    key_positions = {}
    for case_index, angle in enumerate(angles):
        for joint_index, identifier in enumerate(identifiers):
            key_positions[_case_key(case_index, angle, identifier)] = (
                case_index, joint_index
            )
    unknown_keys = [
        value for value in supplied_keys if value not in key_positions
    ]
    if unknown_keys:
        raise ValueError(
            "SlidingCaseKeys contains unknown case keys: {0}".format(
                unknown_keys
            )
        )

    friction = _positive_number(friction_coefficient, "FrictionCoefficient")
    sliding_factor = _positive_number(
        sliding_safety_factor, "SlidingSafetyFactor"
    )
    supplied_key_set = set(supplied_keys)
    sliding_indices_by_case = [set() for angle in angles]
    for key in supplied_key_set:
        case_index, joint_index = key_positions[key]
        sliding_indices_by_case[case_index].add(joint_index)
    prepared_normals = []
    for plane in planes:
        normal = rg.Vector3d(plane.Normal)
        if not normal.Unitize():
            raise ValueError("JointPlanes contains an invalid normal.")
        prepared_normals.append(normal)
    envelopes = {}
    failures = []
    total_deficit_case_rows = 0

    for case_index, angle in enumerate(angles):
        start = case_index * joint_count
        end = start + joint_count
        if row_ids[start:end] != identifiers:
            failures.append((case_index, angle, "case_joint_id_order_mismatch"))
            continue
        case_sliding_ids = [
            identifiers[joint_index]
            for joint_index in sorted(sliding_indices_by_case[case_index])
        ]
        result = create_sliding_deficit(
            planes,
            identifiers,
            resultants[start:end],
            case_sliding_ids,
            friction,
            sliding_factor,
            prepared_normals,
        )
        if not result[3]:
            failures.append((case_index, angle, result[4]))
            continue
        total_deficit_case_rows += len(result[0])
        for identifier, vector, demand in zip(
            result[0], result[1], result[2]
        ):
            joint_index = identifier_indices[identifier]
            plane = planes[joint_index]
            resultant = resultants[start + joint_index]
            normal = prepared_normals[joint_index]
            normal_force = rg.Vector3d.Multiply(resultant, normal)
            tangential_force = (
                resultant - normal * normal_force
            ).Length
            sliding_utilisation = (
                sliding_factor * tangential_force
                / (friction * normal_force)
            )
            candidate = {
                "case_index": case_index,
                "angle": angle,
                "vector": vector,
                "demand": demand,
                "normal": normal_force,
                "tangential": tangential_force,
                "sliding_utilisation": sliding_utilisation,
                "pressure_centroid": centroids[start + joint_index],
            }
            current = envelopes.get(identifier)
            if (
                current is None
                or demand > current["demand"] + DEMAND_TOLERANCE_KN
                or (
                    abs(demand - current["demand"]) <= DEMAND_TOLERANCE_KN
                    and angle < current["angle"]
                )
            ):
                envelopes[identifier] = candidate

    deficit_ids = [
        identifier for identifier in identifiers if identifier in envelopes
    ]
    demand_vectors = [envelopes[value]["vector"] for value in deficit_ids]
    demand_forces = [envelopes[value]["demand"] for value in deficit_ids]
    governing_angles = [envelopes[value]["angle"] for value in deficit_ids]
    governing_normals = [envelopes[value]["normal"] for value in deficit_ids]
    governing_tangentials = [
        envelopes[value]["tangential"] for value in deficit_ids
    ]
    governing_sliding_utilisations = [
        envelopes[value]["sliding_utilisation"] for value in deficit_ids
    ]
    governing_pressure_centroids = [
        envelopes[value]["pressure_centroid"] for value in deficit_ids
    ]
    one_to_one = (
        len(deficit_ids)
        == len(demand_vectors)
        == len(demand_forces)
        == len(governing_angles)
        == len(governing_normals)
        == len(governing_tangentials)
        == len(governing_sliding_utilisations)
        == len(governing_pressure_centroids)
    )
    valid = not failures and one_to_one
    report = (
        f"stage=directional_deficit_demand_envelope; angle_cases={len(angles)}; "
        f"joints_per_case={joint_count}; deficit_case_rows={total_deficit_case_rows}; "
        f"envelope_deficit_joints={len(deficit_ids)}; deficit_joint_ids={deficit_ids}; "
        f"deficit_demand_force_min_max_kN="
        f"({min(demand_forces) if demand_forces else 0.0:.6g},"
        f"{max(demand_forces) if demand_forces else 0.0:.6g}); "
        f"governing_angle_min_max_deg="
        f"({min(governing_angles) if governing_angles else 0.0:.6g},"
        f"{max(governing_angles) if governing_angles else 0.0:.6g}); "
        f"case_failures={failures}; one_to_one={one_to_one}; valid={valid}; "
        "envelope_model=max_per_joint_of_final_deficit_demand_across_angles; "
        "paired_state=angle_normal_tangent_vector_from_same_governing_case; "
        "per_joint_governing_angles_do_not_form_one_simultaneous_global_case; "
        "pressure_is_assessed_separately"
    )
    return (
        deficit_ids,
        demand_vectors,
        demand_forces,
        governing_angles,
        governing_normals,
        governing_tangentials,
        governing_sliding_utilisations,
        governing_pressure_centroids,
        valid,
        report,
    )


try:
    _component_inputs = (
        JointPlanes,
        JointIDs,
        CaseAngles,
        CaseJointIDs,
        ContactResultants,
        PressureCentroids,
        SlidingCaseKeys,
        FrictionCoefficient,
        SlidingSafetyFactor,
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        DeficitJointIDs,
        DeficitVectors,
        SlidingDeficits,
        GoverningHorizontalAngles,
        GoverningNormalForces,
        GoverningTangentialForces,
        GoverningSlidingUtilisations,
        GoverningPressureCentroids,
        Valid,
        Report,
    ) = create_directional_sliding_deficit_envelope(*_component_inputs)
