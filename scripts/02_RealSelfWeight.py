"""Stage 02 - lump the real solid self-weight onto a radial thrust network.

Grasshopper Python inputs:
    WholeSolid    Brep         Closed union from 01_ParametricGeometry
    CenterLines   list[Curve]  N crown-to-support centre lines
    UnitWeight    float        Stone unit weight in kN/m3 (Item, optional)

Use Item Access for WholeSolid and UnitWeight, and List Access for CenterLines.

Outputs:
    Nodes, Edges, NodalLoads, Valid, Report

The solid is sampled only once, after all crown overlaps and fillets have been
unified. Every occupied sample is assigned to its nearest network node, so the
shared crown volume cannot be counted once per arch. The sampled tributary
volumes are scaled to the exact Brep volume before loads are produced.

The expensive exact Brep classification is cached only when a strong geometry
signature can be built. Unchanged Grasshopper solutions reuse numeric tributary
volumes; any geometry, network, tolerance or sampling change invalidates the
cache. A missing signature disables caching rather than using a weak fallback.

UnitWeight defaults to the current research value when its optional input is
empty. The geometry is assumed to use metres, and loads are returned in kN.
"""

import math

import Rhino
import Rhino.Geometry as rg
import numpy as np

try:
    import scriptcontext as sc
except ImportError:
    sc = None


DEFAULT_SEGMENTS_PER_ARCH = 12
DEFAULT_UNIT_WEIGHT = 25.5  # kN/m3; current medium/high-density limestone study.
DEFAULT_TARGET_SAMPLES = 70000
DEFAULT_MIN_GRID_AXIS = 8
DEFAULT_MAX_GRID_AXIS = 72
NEAREST_NODE_CHUNK_SIZE = 2048
CACHE_ENTRY_LIMIT = 4
CACHE_VERSION = "stage02_exact_brep_v2"
CACHE_STICKY_KEY = "kangaroo2_grasshopper_stage02_exact_selfweight_cache_v2"


try:
    _MODULE_CACHE
except NameError:
    _MODULE_CACHE = {}


def _as_list(value):
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return [value]


def _positive_number(value, default, name):
    if value is None:
        return default
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError("{0} must be a finite positive number.".format(name))
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("{0} must be a finite positive number.".format(name))
    return number


def _tolerance():
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return 1.0e-3
    return document.ModelAbsoluteTolerance


def _network_from_curves(curves, segments, tolerance):
    nodes = []
    edges = []
    crown = None

    for arch_index, curve in enumerate(curves):
        if curve is None or not curve.IsValid:
            raise ValueError(
                "CenterLines contains an invalid curve at index {0}."
                .format(arch_index)
            )
        parameters = curve.DivideByCount(segments, True)
        if parameters is None or len(parameters) != segments + 1:
            raise RuntimeError(
                "Could not divide CenterLines[{0}] into {1} segments."
                .format(arch_index, segments)
            )
        points = [curve.PointAt(parameter) for parameter in parameters]

        if crown is None:
            crown = points[0]
            nodes.append(rg.Point3d(crown))
        elif points[0].DistanceTo(crown) > max(
            10.0 * tolerance, 1.0e-7
        ):
            raise ValueError(
                "All CenterLines must start at one shared crown point."
            )

        previous = 0
        for point in points[1:]:
            current = len(nodes)
            nodes.append(rg.Point3d(point))
            edges.append(rg.Line(nodes[previous], nodes[current]))
            previous = current

    return nodes, edges


def _grid_counts(box, target_samples):
    dx = max(box.Max.X - box.Min.X, 1.0e-9)
    dy = max(box.Max.Y - box.Min.Y, 1.0e-9)
    dz = max(box.Max.Z - box.Min.Z, 1.0e-9)
    base = (float(target_samples) / (dx * dy * dz)) ** (1.0 / 3.0)

    counts = []
    for length in (dx, dy, dz):
        count = int(round(length * base))
        count = max(DEFAULT_MIN_GRID_AXIS, count)
        count = min(DEFAULT_MAX_GRID_AXIS, count)
        counts.append(count)
    return counts


def _nearest_node_index(point, nodes):
    nearest_index = 0
    nearest_distance = float("inf")
    for index, node in enumerate(nodes):
        dx = point.X - node.X
        dy = point.Y - node.Y
        dz = point.Z - node.Z
        distance = dx * dx + dy * dy + dz * dz
        if distance < nearest_distance:
            nearest_distance = distance
            nearest_index = index
    return nearest_index


def _geometry_crc(geometry):
    try:
        return int(geometry.DataCRC(0))
    except Exception:
        return None


def _cache_key(solid, nodes, target_samples, tolerance):
    crc = _geometry_crc(solid)
    if crc is None:
        return None
    box = solid.GetBoundingBox(True)
    try:
        topology = (
            int(solid.Faces.Count),
            int(solid.Edges.Count),
            int(solid.Vertices.Count),
        )
    except Exception:
        return None
    volume = solid.GetVolume()
    if volume <= 0.0:
        return None
    return (
        CACHE_VERSION,
        crc,
        topology,
        (
            float(box.Min.X), float(box.Min.Y), float(box.Min.Z),
            float(box.Max.X), float(box.Max.Y), float(box.Max.Z),
        ),
        float(volume),
        tuple(
            (float(node.X), float(node.Y), float(node.Z))
            for node in nodes
        ),
        int(target_samples),
        float(tolerance),
        int(DEFAULT_MIN_GRID_AXIS),
        int(DEFAULT_MAX_GRID_AXIS),
    )


def _cache_store():
    if sc is not None:
        return sc.sticky
    return _MODULE_CACHE


def _cached_tributary_result(key):
    if key is None:
        return None
    store = _cache_store()
    bucket = store.get(CACHE_STICKY_KEY)
    if not isinstance(bucket, dict):
        return None
    entries = bucket.get("entries")
    if not isinstance(entries, list):
        return None
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict) or entry.get("key") != key:
            continue
        payload = entry.get("payload")
        if not isinstance(payload, tuple) or len(payload) != 5:
            return None
        if index:
            entries.insert(0, entries.pop(index))
        volumes, grid, inside_count, exact_volume, sampled_volume = payload
        return (
            list(volumes),
            tuple(grid),
            int(inside_count),
            float(exact_volume),
            float(sampled_volume),
        )
    return None


def _store_tributary_result(key, result):
    if key is None:
        return
    volumes, grid, inside_count, exact_volume, sampled_volume = result
    store = _cache_store()
    bucket = store.get(CACHE_STICKY_KEY)
    if not isinstance(bucket, dict) or not isinstance(
        bucket.get("entries"), list
    ):
        bucket = {"entries": []}
    entries = [
        entry for entry in bucket["entries"]
        if isinstance(entry, dict) and entry.get("key") != key
    ]
    entry = {
        "key": key,
        "payload": (
            tuple(float(value) for value in volumes),
            tuple(int(value) for value in grid),
            int(inside_count),
            float(exact_volume),
            float(sampled_volume),
        ),
    }
    entries.insert(0, entry)
    bucket["entries"] = entries[:CACHE_ENTRY_LIMIT]
    store[CACHE_STICKY_KEY] = bucket


def _nearest_node_indices(coordinates, nodes):
    node_coordinates = np.asarray(
        [(node.X, node.Y, node.Z) for node in nodes], dtype=float
    )
    sample_coordinates = np.asarray(coordinates, dtype=float)
    nearest = []
    for start in range(0, len(sample_coordinates), NEAREST_NODE_CHUNK_SIZE):
        chunk = sample_coordinates[start:start + NEAREST_NODE_CHUNK_SIZE]
        delta_x = chunk[:, None, 0] - node_coordinates[None, :, 0]
        delta_y = chunk[:, None, 1] - node_coordinates[None, :, 1]
        delta_z = chunk[:, None, 2] - node_coordinates[None, :, 2]
        distances = (
            delta_x * delta_x + delta_y * delta_y + delta_z * delta_z
        )
        nearest.extend(np.argmin(distances, axis=1).tolist())
    return nearest


def _sample_tributary_volumes(solid, nodes, target_samples, tolerance):
    key = _cache_key(solid, nodes, target_samples, tolerance)
    cached = _cached_tributary_result(key)
    if cached is not None:
        return cached + (True, key is not None)

    box = solid.GetBoundingBox(True)
    nx, ny, nz = _grid_counts(box, target_samples)
    dx = (box.Max.X - box.Min.X) / float(nx)
    dy = (box.Max.Y - box.Min.Y) / float(ny)
    dz = (box.Max.Z - box.Min.Z) / float(nz)
    voxel_volume = dx * dy * dz

    inside_coordinates = []
    for ix in range(nx):
        x = box.Min.X + (ix + 0.5) * dx
        for iy in range(ny):
            y = box.Min.Y + (iy + 0.5) * dy
            for iz in range(nz):
                z = box.Min.Z + (iz + 0.5) * dz
                point = rg.Point3d(x, y, z)
                if not solid.IsPointInside(point, tolerance, False):
                    continue
                inside_coordinates.append((x, y, z))

    inside_count = len(inside_coordinates)

    if inside_count == 0:
        raise RuntimeError(
            "No interior volume samples were found. Check model units and "
            "whether WholeSolid is a closed solid."
        )

    volumes = [0.0 for node in nodes]
    for node_index in _nearest_node_indices(inside_coordinates, nodes):
        volumes[node_index] += voxel_volume

    exact_volume = solid.GetVolume()
    sampled_volume = sum(volumes)
    if exact_volume <= 0.0 or sampled_volume <= 0.0:
        raise RuntimeError("WholeSolid has no positive volume.")
    scale = exact_volume / sampled_volume
    volumes = [volume * scale for volume in volumes]
    result = (
        volumes, (nx, ny, nz), inside_count, exact_volume, sampled_volume
    )
    _store_tributary_result(key, result)
    return result + (False, key is not None)


def create_real_selfweight(whole_solid, center_lines, unit_weight=None):
    curves = _as_list(center_lines)
    if whole_solid is None:
        raise ValueError("WholeSolid is required.")
    if not whole_solid.IsValid or not whole_solid.IsSolid:
        raise ValueError("WholeSolid must be one valid closed Brep.")
    if len(curves) < 3:
        raise ValueError("CenterLines must contain at least three curves.")

    weight_density = _positive_number(
        unit_weight, DEFAULT_UNIT_WEIGHT, "UnitWeight"
    )
    segments = DEFAULT_SEGMENTS_PER_ARCH
    tolerance = _tolerance()
    nodes, edges = _network_from_curves(
        curves, segments, tolerance
    )
    (
        tributary_volumes,
        grid,
        inside_count,
        exact_volume,
        sampled_volume,
        cache_hit,
        cache_available,
    ) = _sample_tributary_volumes(
        whole_solid,
        nodes,
        DEFAULT_TARGET_SAMPLES,
        tolerance,
    )

    nodal_weights = [
        weight_density * volume for volume in tributary_volumes
    ]
    nodal_loads = [
        rg.Vector3d(0.0, 0.0, -weight) for weight in nodal_weights
    ]
    total_weight = sum(nodal_weights)
    expected_weight = weight_density * exact_volume
    weight_error = abs(total_weight - expected_weight)

    crown_weight = nodal_weights[0]
    valid = (
        len(nodes) == 1 + len(curves) * segments
        and len(edges) == len(curves) * segments
        and len(nodal_loads) == len(nodes)
        and all(weight >= 0.0 for weight in nodal_weights)
        and weight_error <= max(1.0e-9, expected_weight * 1.0e-10)
    )

    report = (
        "arches={0}; segments_per_arch={1}; nodes={2}; edges={3}; "
        "grid={4}x{5}x{6}; inside_samples={7}; "
        "exact_volume={8:.6g}m3; sampled_volume={9:.6g}m3; "
        "unit_weight={10:.6g}kN/m3; total_weight={11:.6g}kN; "
        "crown_weight_once={12:.6g}kN; "
        "weight_error={13:.6g}kN; cache_available={14}; cache_hit={15}; "
        "cache_entry_limit={16}; "
        "nearest_assignment=chunked_numpy_first_argmin; valid={17}"
    ).format(
        len(curves),
        segments,
        len(nodes),
        len(edges),
        grid[0],
        grid[1],
        grid[2],
        inside_count,
        exact_volume,
        sampled_volume,
        weight_density,
        total_weight,
        crown_weight,
        weight_error,
        cache_available,
        cache_hit,
        CACHE_ENTRY_LIMIT,
        valid,
    )
    return nodes, edges, nodal_loads, valid, report


try:
    _component_inputs = (WholeSolid, CenterLines, UnitWeight)
except NameError:
    try:
        _component_inputs = (WholeSolid, CenterLines, None)
    except NameError:
        _component_inputs = None

if _component_inputs is not None:
    Nodes, Edges, NodalLoads, Valid, Report = create_real_selfweight(
        *_component_inputs
    )
