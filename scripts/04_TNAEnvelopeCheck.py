#! python 3
"""Stage 04 - compression-only TNA equilibrium and solid containment check.

Grasshopper Python inputs:
    WholeSolid    Brep          From 01_ParametricGeometry
    InitialEdges  list[Line]    ThrustEdges from 03_KangarooInitialForm
    NodalLoads    list[Vector3d] From 02_RealSelfWeight

Use Item Access for WholeSolid and List Access for the other two inputs.

Outputs:
    TNAEdges, ThrustLines, OutsidePoints, Valid, Report, Advice

COMPAS TNA's force-density formulation is adapted to the project's radial
tree topology. Standard FormDiagram objects require faces, whereas this
structure is N one-dimensional half arches meeting at one node. Horizontal
compression magnitudes are optimised subject to exact crown equilibrium.
Vertical equilibrium is then solved for the real nonuniform nodal loads.
The crown is a free internal equilibrium node; only the outer leaves are
supports. The resulting network is sampled against the actual Rhino Brep.
"""

import math

import compas
import compas_tna
import numpy as np
import Rhino
import Rhino.Geometry as rg
from compas.matrices import connectivity_matrix
from scipy.optimize import minimize
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


DEFAULT_EDGE_SAMPLES = 12
DEFAULT_EQUILIBRIUM_TOLERANCE = 1.0e-6
DEFAULT_COMPRESSION_EPSILON = 1.0e-8
DEFAULT_MIN_HORIZONTAL_WEIGHT = 1.0e-4
DEFAULT_MAX_HORIZONTAL_WEIGHT = 100.0


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
    indexed_edges = []
    seen_edges = set()

    def vertex_index(point):
        key = _point_key(point, tolerance)
        if key not in key_to_index:
            key_to_index[key] = len(points)
            points.append(rg.Point3d(point))
        return key_to_index[key]

    for line_index, line in enumerate(lines):
        if not line.IsValid or line.Length <= tolerance:
            raise ValueError(
                "InitialEdges contains an invalid line at index {0}."
                .format(line_index)
            )
        a = vertex_index(line.From)
        b = vertex_index(line.To)
        key = (min(a, b), max(a, b))
        if key in seen_edges:
            raise ValueError("InitialEdges contains a duplicate edge.")
        seen_edges.add(key)
        indexed_edges.append(key)
    return points, indexed_edges


def _radial_topology(points, edges):
    adjacency = [[] for _ in points]
    for edge_index, (a, b) in enumerate(edges):
        adjacency[a].append((b, edge_index))
        adjacency[b].append((a, edge_index))

    supports = [
        index for index, neighbours in enumerate(adjacency)
        if len(neighbours) == 1
    ]
    crowns = [
        index for index, neighbours in enumerate(adjacency)
        if len(neighbours) > 2
    ]
    if len(crowns) != 1:
        raise ValueError("Exactly one shared crown node is required.")
    if len(supports) < 3:
        raise ValueError("At least three support leaves are required.")
    if len(edges) != len(points) - 1:
        raise ValueError("The TNA network must be one radial tree.")

    crown = crowns[0]
    branches = []
    branch_edges = []
    for first, first_edge in adjacency[crown]:
        nodes = [crown, first]
        edge_indices = [first_edge]
        previous = crown
        current = first
        while current not in supports:
            following = [
                item for item in adjacency[current]
                if item[0] != previous
            ]
            if len(following) != 1:
                raise ValueError(
                    "Every crown-to-support path must be a simple chain."
                )
            next_node, next_edge = following[0]
            nodes.append(next_node)
            edge_indices.append(next_edge)
            previous, current = current, next_node
        branches.append(nodes)
        branch_edges.append(edge_indices)
    return crown, supports, branches, branch_edges


def _branch_directions(points, crown, branches):
    directions = []
    for branch in branches:
        point = points[branch[1]]
        dx = point.X - points[crown].X
        dy = point.Y - points[crown].Y
        length = math.sqrt(dx * dx + dy * dy)
        if length <= 1.0e-9:
            raise ValueError("A crown edge has zero plan length.")
        directions.append((dx / length, dy / length))
    return np.asarray(directions, dtype=float)


def _balanced_start(directions):
    count = len(directions)
    constraints = [
        {
            "type": "eq",
            "fun": lambda h: float(np.dot(h, directions[:, 0])),
        },
        {
            "type": "eq",
            "fun": lambda h: float(np.dot(h, directions[:, 1])),
        },
        {
            "type": "eq",
            "fun": lambda h: float(np.mean(h) - 1.0),
        },
    ]
    result = minimize(
        lambda h: float(np.sum((h - 1.0) ** 2)),
        np.ones(count),
        method="SLSQP",
        bounds=[
            (
                DEFAULT_MIN_HORIZONTAL_WEIGHT,
                DEFAULT_MAX_HORIZONTAL_WEIGHT,
            )
            for _ in range(count)
        ],
        constraints=constraints,
        options={"maxiter": 500, "ftol": 1.0e-12, "disp": False},
    )
    if not result.success:
        raise RuntimeError(
            "No positive horizontal crown equilibrium was found: {0}"
            .format(result.message)
        )
    return np.asarray(result.x, dtype=float), constraints


def _edge_force_densities(points, edges, branch_edges, weights):
    q = np.zeros(len(edges), dtype=float)
    for branch_index, edge_indices in enumerate(branch_edges):
        horizontal_force = weights[branch_index]
        for edge_index in edge_indices:
            a, b = edges[edge_index]
            dx = points[b].X - points[a].X
            dy = points[b].Y - points[a].Y
            plan_length = math.sqrt(dx * dx + dy * dy)
            if plan_length <= 1.0e-9:
                raise ValueError("A TNA edge has zero plan length.")
            q[edge_index] = -horizontal_force / plan_length
    return q


def _solve_vertical(points, edges, supports, crown, loads,
                    branch_edges, weights, target_crown):
    vertex_count = len(points)
    fixed = sorted(supports)
    fixed_set = set(fixed)
    free = [
        index for index in range(vertex_count)
        if index not in fixed_set
    ]
    C = connectivity_matrix(edges, rtype="csr")
    q0 = _edge_force_densities(
        points, edges, branch_edges, weights
    )
    Q0 = diags(q0)
    K0 = C.transpose().dot(Q0).dot(C)
    pz = np.asarray([load.Z for load in loads], dtype=float)
    Kff0 = K0[free, :][:, free]
    Kfc0 = K0[free, :][:, fixed]
    fixed_z = np.zeros(len(fixed), dtype=float)
    z0 = np.zeros(vertex_count, dtype=float)
    z0[free] = spsolve(
        Kff0, pz[free] - Kfc0.dot(fixed_z)
    )

    if z0[crown] <= 1.0e-9 or target_crown <= 1.0e-9:
        raise RuntimeError(
            "The compression solution produced no positive crown height."
        )
    scale = z0[crown] / target_crown
    q = q0 * scale
    Q = diags(q)
    K = C.transpose().dot(Q).dot(C)
    Kff = K[free, :][:, free]
    Kfc = K[free, :][:, fixed]
    z = np.zeros(vertex_count, dtype=float)
    z[free] = spsolve(Kff, pz[free] - Kfc.dot(fixed_z))
    return z, q, K, free, fixed, scale


def _optimise_shape(points, edges, supports, crown, branch_edges,
                    loads, initial_z, directions):
    start, constraints = _balanced_start(directions)
    target_crown = initial_z[crown]
    height_scale = max(target_crown, 1.0)

    def objective(weights):
        try:
            z = _solve_vertical(
                points,
                edges,
                supports,
                crown,
                loads,
                branch_edges,
                weights,
                target_crown,
            )[0]
        except Exception:
            return 1.0e12
        shape_error = np.mean(
            ((z - initial_z) / height_scale) ** 2
        )
        regularisation = 1.0e-4 * np.mean((weights - start) ** 2)
        return float(shape_error + regularisation)

    result = minimize(
        objective,
        start,
        method="SLSQP",
        bounds=[
            (
                DEFAULT_MIN_HORIZONTAL_WEIGHT,
                DEFAULT_MAX_HORIZONTAL_WEIGHT,
            )
            for _ in range(len(start))
        ],
        constraints=constraints,
        options={"maxiter": 600, "ftol": 1.0e-12, "disp": False},
    )
    weights = (
        np.asarray(result.x, dtype=float)
        if result.success
        else start
    )
    solution = _solve_vertical(
        points,
        edges,
        supports,
        crown,
        loads,
        branch_edges,
        weights,
        target_crown,
    )
    return weights, result, solution


def _containment(solid, points, edges, tolerance):
    outside = []
    outside_edge_indices = []
    for edge_index, (a, b) in enumerate(edges):
        start = points[a]
        end = points[b]
        for sample_index in range(DEFAULT_EDGE_SAMPLES + 1):
            fraction = (
                float(sample_index) / float(DEFAULT_EDGE_SAMPLES)
            )
            point = rg.Point3d(
                start.X + (end.X - start.X) * fraction,
                start.Y + (end.Y - start.Y) * fraction,
                start.Z + (end.Z - start.Z) * fraction,
            )
            if not solid.IsPointInside(point, tolerance, False):
                outside.append(point)
                outside_edge_indices.append(edge_index)
    sample_count = len(edges) * (DEFAULT_EDGE_SAMPLES + 1)
    return outside, outside_edge_indices, sample_count


def _make_baseline_advice(valid, compression_only, equilibrium,
                          outside_edge_indices, branch_edges):
    """Explain selected-network equilibrium and sampled containment only."""
    if valid:
        return (
            "SELECTED-NETWORK PASS: compression-only equilibrium is satisfied; "
            "the sampled points of every discrete TNA edge lie inside WholeSolid. "
            "This is sampled evidence, not continuous containment or structural safety. "
            "Actual joint pressure and sliding deficit require separate assessment."
        )

    messages = []
    if not equilibrium:
        messages.append(
            "STOP: the current network does not satisfy free-node "
            "equilibrium. Do not use it for segmentation. Revise the "
            "support layout, RotationAngles or network topology and rerun "
            "Stages 01-04."
        )
    if not compression_only:
        messages.append(
            "SELECTED-NETWORK COMPRESSION FAILURE: the current network "
            "contains tensile force. Omit downstream assessment for this candidate "
            "and return to designer revision. This does not rule out other networks."
        )
    if outside_edge_indices:
        edge_to_branch = {}
        for branch_index, edge_indices in enumerate(branch_edges):
            for edge_index in edge_indices:
                edge_to_branch[edge_index] = branch_index
        branch_counts = {}
        for edge_index in outside_edge_indices:
            branch_index = edge_to_branch.get(edge_index, -1)
            branch_counts[branch_index] = (
                branch_counts.get(branch_index, 0) + 1
            )
        affected = [
            branch_index + 1
            for branch_index, count in sorted(
                branch_counts.items(),
                key=lambda item: item[1],
                reverse=True,
            )
            if branch_index >= 0
        ]
        messages.append(
            "BASELINE CONTAINMENT WARNING: {0} sampled points leave "
            "WholeSolid; affected arch indices={1}. The equilibrium "
            "solution exists, but sampled containment of this selected network "
            "fails. Omit downstream results; the designer may revise H, Spans "
            "or CrownDepth and regenerate the complete candidate. This "
            "does not prove that no other admissible network exists.".format(
                len(outside_edge_indices), affected
            )
        )
    return " ".join(messages)


def _branch_polylines(points, branches):
    return [
        rg.PolylineCurve([points[index] for index in branch])
        for branch in branches
    ]


def create_tna_envelope_check(whole_solid, initial_edges, nodal_loads):
    lines = _as_list(initial_edges)
    loads = [rg.Vector3d(load) for load in _as_list(nodal_loads)]
    if whole_solid is None:
        raise ValueError("WholeSolid is required.")
    if not whole_solid.IsValid or not whole_solid.IsSolid:
        raise ValueError("WholeSolid must be one valid closed Brep.")
    if len(lines) < 3:
        raise ValueError("InitialEdges must contain a radial network.")

    tolerance = _document_tolerance()
    initial_points, indexed_edges = _ordered_network(
        lines, max(tolerance * 0.1, 1.0e-7)
    )
    if len(loads) != len(initial_points):
        raise ValueError(
            "NodalLoads count must match the unique InitialEdges nodes."
        )
    (
        crown,
        supports,
        branches,
        branch_edges,
    ) = _radial_topology(initial_points, indexed_edges)
    directions = _branch_directions(
        initial_points, crown, branches
    )
    initial_z = np.asarray(
        [point.Z for point in initial_points], dtype=float
    )
    weights, optimisation, solution = _optimise_shape(
        initial_points,
        indexed_edges,
        supports,
        crown,
        branch_edges,
        loads,
        initial_z,
        directions,
    )
    z, q, K, free, fixed, force_scale = solution

    xyz = np.asarray(
        [[point.X, point.Y, z[index]]
         for index, point in enumerate(initial_points)],
        dtype=float,
    )
    p = np.zeros((len(initial_points), 3), dtype=float)
    p[:, 2] = np.asarray([load.Z for load in loads], dtype=float)
    residuals = K.dot(xyz) - p
    free_residual = float(np.linalg.norm(residuals[free, :]))
    crown_residual = float(np.linalg.norm(residuals[crown, :]))
    tna_points = [
        rg.Point3d(row[0], row[1], row[2]) for row in xyz
    ]
    tna_edges = [
        rg.Line(tna_points[a], tna_points[b])
        for a, b in indexed_edges
    ]
    thrust_lines = _branch_polylines(tna_points, branches)
    (
        outside_points,
        outside_edge_indices,
        containment_samples,
    ) = _containment(
        whole_solid, tna_points, indexed_edges, tolerance
    )

    edge_forces = []
    for edge_index, (a, b) in enumerate(indexed_edges):
        length = tna_points[a].DistanceTo(tna_points[b])
        edge_forces.append(-q[edge_index] * length)
    compression_only = (
        all(value > DEFAULT_COMPRESSION_EPSILON for value in edge_forces)
        and all(value < -DEFAULT_COMPRESSION_EPSILON for value in q)
    )
    horizontal_balance = math.sqrt(
        float(np.dot(weights, directions[:, 0])) ** 2
        + float(np.dot(weights, directions[:, 1])) ** 2
    )
    inside = len(outside_points) == 0
    total_weight = -sum(load.Z for load in loads)
    residual_tolerance = max(
        DEFAULT_EQUILIBRIUM_TOLERANCE,
        total_weight * 1.0e-6,
    )
    equilibrium = (
        free_residual <= residual_tolerance
        and crown_residual <= residual_tolerance
        and horizontal_balance <= DEFAULT_EQUILIBRIUM_TOLERANCE
    )
    valid = compression_only and equilibrium and inside
    advice = _make_baseline_advice(
        valid,
        compression_only,
        equilibrium,
        outside_edge_indices,
        branch_edges,
    )
    reactions = [
        (
            round(float(residuals[index, 0]), 6),
            round(float(residuals[index, 1]), 6),
            round(float(residuals[index, 2]), 6),
        )
        for index in fixed
    ]
    report = (
        "engine=COMPAS_TNA_{0}_radial_graph_adapter; "
        "compas={1}; nodes={2}; edges={3}; arches={4}; "
        "supports={5}; crown={6}; real_nodal_loads=True; "
        "total_weight={7:.6g}kN; optimisation_success={8}; "
        "optimisation_iterations={9}; horizontal_weights={10}; "
        "force_scale={11:.6g}; min_compression={12:.6g}kN; "
        "max_compression={13:.6g}kN; compression_only={14}; "
        "horizontal_balance={15:.6g}; free_residual={16:.6g}; "
        "crown_residual={17:.6g}; centre_is_free=True; "
        "containment_samples={18}; outside_samples={19}; "
        "inside_whole_solid={20}; support_reactions_kN={21}; "
        "valid={22}"
    ).format(
        compas_tna.__version__,
        compas.__version__,
        len(tna_points),
        len(tna_edges),
        len(branches),
        supports,
        crown,
        total_weight,
        bool(optimisation.success),
        int(getattr(optimisation, "nit", 0)),
        [round(float(value), 6) for value in weights],
        force_scale,
        min(edge_forces),
        max(edge_forces),
        compression_only,
        horizontal_balance,
        free_residual,
        crown_residual,
        containment_samples,
        len(outside_points),
        inside,
        reactions,
        valid,
    )
    return (
        tna_edges,
        thrust_lines,
        outside_points,
        valid,
        report,
        advice,
    )


try:
    _component_inputs = (WholeSolid, InitialEdges, NodalLoads)
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    TNAEdges, ThrustLines, OutsidePoints, Valid, Report, Advice = (
        create_tna_envelope_check(*_component_inputs)
    )
