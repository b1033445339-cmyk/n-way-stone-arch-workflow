#! python 3
"""Stage 07A - exact block self-weight and actual joint adjacency.

Grasshopper Python inputs:
    Blocks          list[Brep]   Blocks from Stage 06 (List)
    JointPlanes     list[Plane]  JointPlanes from Stage 06 (List)
    JointProfiles   list[Curve]  JointProfiles from Stage 06 (List)
    BlockIDs        list[str]    BlockIDs from Stage 06 (List)
    UnitWeight      float        Stone unit weight in kN/m3 (Item, optional)

Outputs:
    BlockCentroids, BlockLoads, JointIDs, JointBlockPairs, Valid, Report

This first part of Stage 07 establishes exact block weights and the topology of
the already segmented assembly. It does not check contact pressure, opening,
friction, crushing or dowel capacity.
"""

import math
import re

import Rhino
import Rhino.Geometry as rg


DEFAULT_UNIT_WEIGHT = 25.5  # kN/m3; match the current Stage 02 limestone study.


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


def _coerce_brep(value):
    if isinstance(value, rg.Brep):
        return value
    if hasattr(value, "Value") and isinstance(value.Value, rg.Brep):
        return value.Value
    return None


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


def _positive_number(value, default, name):
    if value is None:
        return default
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("{0} must be a finite positive number.".format(name))
    return number


def _volume_centroid(block):
    properties = rg.VolumeMassProperties.Compute(block)
    if properties is None:
        raise RuntimeError("Could not compute a block volume and centroid.")
    try:
        return abs(properties.Volume), rg.Point3d(properties.Centroid)
    finally:
        properties.Dispose()


def _profile_centroid(profile):
    properties = rg.AreaMassProperties.Compute(profile)
    if properties is None:
        raise RuntimeError("Could not compute a joint-profile centroid.")
    try:
        return rg.Point3d(properties.Centroid), abs(properties.Area)
    finally:
        properties.Dispose()


def _containing_blocks(blocks, point, tolerance):
    return [
        index for index, block in enumerate(blocks)
        if block.IsPointInside(point, tolerance, False)
    ]


def _joint_pair(blocks, profile, plane, tolerance):
    centroid, area = _profile_centroid(profile)
    normal = rg.Vector3d(plane.Normal)
    if not normal.Unitize():
        raise ValueError("A JointPlane has an invalid normal.")

    base_offset = max(20.0 * tolerance, 1.0e-4)
    attempts = []
    for factor in (1.0, 2.0, 5.0, 10.0):
        offset = base_offset * factor
        plus = _containing_blocks(
            blocks, centroid + normal * offset, tolerance
        )
        minus = _containing_blocks(
            blocks, centroid - normal * offset, tolerance
        )
        attempts.append((offset, plus, minus))
        if len(plus) == 1 and len(minus) == 1 and plus[0] != minus[0]:
            return plus[0], minus[0], area, offset, attempts
    return None, None, area, None, attempts


def _parse_arch_block(identifier):
    match = re.fullmatch(r"arch\[(\d+)\]/block\[(\d+)\]", identifier)
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def _joint_id(first_id, second_id):
    first = _parse_arch_block(first_id)
    second = _parse_arch_block(second_id)
    first_is_center = first_id.startswith("center/")
    second_is_center = second_id.startswith("center/")

    if first_is_center and second is not None and second[1] == 1:
        return "arch[{0}]/joint[1]".format(second[0])
    if second_is_center and first is not None and first[1] == 1:
        return "arch[{0}]/joint[1]".format(first[0])
    if first is None or second is None or first[0] != second[0]:
        return None
    if abs(first[1] - second[1]) != 1:
        return None
    return "arch[{0}]/joint[{1}]".format(
        first[0], max(first[1], second[1])
    )


def _connected_block_count(block_count, pairs):
    if block_count == 0:
        return 0
    adjacency = [[] for _ in range(block_count)]
    for first, second in pairs:
        adjacency[first].append(second)
        adjacency[second].append(first)
    visited = {0}
    pending = [0]
    while pending:
        current = pending.pop()
        for neighbour in adjacency[current]:
            if neighbour not in visited:
                visited.add(neighbour)
                pending.append(neighbour)
    return len(visited)


def create_block_weight_and_joint_adjacency(
    blocks,
    joint_planes,
    joint_profiles,
    block_ids,
    unit_weight=None,
):
    block_values = []
    for index, value in enumerate(_as_list(blocks)):
        block = _coerce_brep(value)
        if block is None or not block.IsValid or not block.IsSolid:
            raise ValueError(
                "Blocks contains an invalid closed Brep at index {0}."
                .format(index)
            )
        block_values.append(block)
    if not block_values:
        raise ValueError("Blocks must contain valid closed Breps.")

    identifiers = [str(value) for value in _as_list(block_ids)]
    if len(identifiers) != len(block_values):
        raise ValueError("BlockIDs must correspond one-to-one with Blocks.")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("BlockIDs must be unique.")
    center_ids = [value for value in identifiers if value.startswith("center/")]
    if len(center_ids) != 1:
        raise ValueError("BlockIDs must contain exactly one center block.")

    planes = []
    for index, value in enumerate(_as_list(joint_planes)):
        plane = _coerce_plane(value)
        if plane is None or not plane.IsValid:
            raise ValueError(
                "JointPlanes contains an invalid plane at index {0}."
                .format(index)
            )
        planes.append(plane)

    profiles = []
    for index, value in enumerate(_as_list(joint_profiles)):
        profile = _coerce_curve(value)
        if profile is None or not profile.IsValid or not profile.IsClosed:
            raise ValueError(
                "JointProfiles contains an invalid closed curve at index {0}."
                .format(index)
            )
        profiles.append(profile)
    if len(planes) != len(profiles):
        raise ValueError(
            "JointPlanes and JointProfiles must correspond one-to-one."
        )

    weight_density = _positive_number(
        unit_weight, DEFAULT_UNIT_WEIGHT, "UnitWeight"
    )
    tolerance = _document_tolerance()
    centroids = []
    loads = []
    volumes = []
    for block in block_values:
        volume, centroid = _volume_centroid(block)
        if volume <= 0.0:
            raise RuntimeError("A block has non-positive volume.")
        volumes.append(volume)
        centroids.append(centroid)
        loads.append(rg.Vector3d(0.0, 0.0, -weight_density * volume))

    joint_ids = []
    pair_indices = []
    pair_ids = []
    areas = []
    offsets = []
    failures = []
    for index, (plane, profile) in enumerate(zip(planes, profiles)):
        first, second, area, offset, attempts = _joint_pair(
            block_values, profile, plane, tolerance
        )
        areas.append(area)
        offsets.append(offset)
        if first is None or second is None:
            failures.append((index, attempts))
            continue
        first_id = identifiers[first]
        second_id = identifiers[second]
        identifier = _joint_id(first_id, second_id)
        if identifier is None:
            failures.append((index, first_id, second_id, "non_radial_pair"))
            continue
        joint_ids.append(identifier)
        pair_indices.append((first, second))
        pair_ids.append((first_id, second_id))

    expected_joint_count = len(block_values) - 1
    unique_pairs = {tuple(sorted(pair)) for pair in pair_indices}
    connected_count = _connected_block_count(len(block_values), pair_indices)
    topology_valid = (
        not failures
        and len(planes) == expected_joint_count
        and len(pair_indices) == len(planes)
        and len(unique_pairs) == len(pair_indices)
        and len(set(joint_ids)) == len(joint_ids)
        and connected_count == len(block_values)
    )
    total_volume = sum(volumes)
    total_weight = sum(-load.Z for load in loads)
    expected_weight = weight_density * total_volume
    weight_error = abs(total_weight - expected_weight)
    weight_tolerance = max(expected_weight * 1.0e-12, 1.0e-9)
    valid = topology_valid and weight_error <= weight_tolerance
    report = (
        f"stage=block_weight_and_joint_adjacency; blocks={len(block_values)}; "
        f"joints={len(planes)}; unit_weight_kN_m3={weight_density:.6g}; "
        f"total_volume_m3={total_volume:.6g}; "
        f"total_weight_kN={total_weight:.6g}; "
        f"weight_balance_error_kN={weight_error:.6g}; "
        f"joint_pairs_resolved={len(pair_indices)}/{len(planes)}; "
        f"unique_joint_pairs={len(unique_pairs)}; "
        f"connected_blocks={connected_count}/{len(block_values)}; "
        f"joint_area_min_max_m2="
        f"({min(areas) if areas else 0.0:.6g},"
        f"{max(areas) if areas else 0.0:.6g}); "
        f"probe_offsets={offsets}; failures={failures}; "
        f"topology_valid={topology_valid}; valid={valid}; "
        "scope=weights_and_adjacency_only; "
        "contact_friction_crushing_and_pin_capacity_not_checked"
    )
    return centroids, loads, joint_ids, pair_ids, valid, report


try:
    _component_inputs = (
        Blocks,
        JointPlanes,
        JointProfiles,
        BlockIDs,
        UnitWeight,
    )
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    (
        BlockCentroids,
        BlockLoads,
        JointIDs,
        JointBlockPairs,
        Valid,
        Report,
    ) = create_block_weight_and_joint_adjacency(*_component_inputs)
