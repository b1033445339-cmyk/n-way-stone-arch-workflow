"""Stage 03 - Kangaroo initial thrust form under real nonuniform self-weight.

Grasshopper Python inputs:
    Nodes         list[Point3d]  From 02_RealSelfWeight
    Edges         list[Line]     From 02_RealSelfWeight
    NodalLoads    list[Vector3d] From 02_RealSelfWeight

Use List Access for all three inputs.

Outputs:
    ThrustEdges, ThrustLines, Converged, Report

This is an initial form-finding result, not the final masonry proof. It uses
the actual nonuniform nodal loads and one coupled crown particle. The outer
leaf nodes are supports. The hanging network is solved in Kangaroo, inverted,
and scaled so the shared crown reaches the crown-centre height supplied by
Nodes. COMPAS TNA/TNO will perform the later compression-only and envelope
checks.
"""

import clr

clr.AddReference("KangarooSolver")

import Rhino.Geometry as rg
from System.Collections.Generic import List
from KangarooSolver import IGoal, PhysicalSystem
from KangarooSolver.Goals import Anchor, Spring, Unary


DEFAULT_REST_LENGTH_FACTOR = 0.94
DEFAULT_STIFFNESS_RATIO = 200.0
DEFAULT_MAX_ITERATIONS = 5000
DEFAULT_MAX_STIFFNESS_RUNS = 12
DEFAULT_SOLVER_TOLERANCE = 1.0e-10


def _as_list(value):
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return [value]


def _nearest_index(point, points, tolerance):
    best = -1
    best_distance = float("inf")
    for index, candidate in enumerate(points):
        distance = point.DistanceTo(candidate)
        if distance < best_distance:
            best = index
            best_distance = distance
    if best_distance > tolerance:
        raise ValueError(
            "An Edges endpoint does not match a Nodes point; "
            "distance={0:.6g}.".format(best_distance)
        )
    return best


def _topology(points, edge_lines):
    tolerance = 1.0e-7
    indexed_edges = []
    adjacency = [[] for _ in points]
    seen = set()
    for edge_index, line in enumerate(edge_lines):
        if not line.IsValid or line.Length <= tolerance:
            raise ValueError(
                "Edges contains an invalid line at index {0}."
                .format(edge_index)
            )
        a = _nearest_index(line.From, points, tolerance)
        b = _nearest_index(line.To, points, tolerance)
        if a == b:
            raise ValueError("Edges contains a zero-topology edge.")
        key = (min(a, b), max(a, b))
        if key in seen:
            raise ValueError("Edges contains a duplicate edge.")
        seen.add(key)
        indexed_edges.append(key)
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
        raise ValueError(
            "The network must contain exactly one shared crown node."
        )
    crown = crowns[0]
    if len(supports) < 3:
        raise ValueError("The network must contain at least three supports.")
    if len(indexed_edges) != len(points) - 1:
        raise ValueError("The radial network must be one connected tree.")

    branches = []
    for first, _ in adjacency[crown]:
        branch = [crown, first]
        previous = crown
        current = first
        while current not in supports:
            following = [
                item for item in adjacency[current]
                if item[0] != previous
            ]
            if len(following) != 1:
                raise ValueError(
                    "Each crown-to-support branch must be a simple chain."
                )
            next_node, _ = following[0]
            branch.append(next_node)
            previous, current = current, next_node
        branches.append(branch)

    if len(branches) != len(supports):
        raise ValueError("Every crown branch must end at one support.")
    return indexed_edges, crown, supports, branches


def _initial_hanging_points(source_points, crown, branches):
    target_rise = source_points[crown].Z
    result = [
        rg.Point3d(point.X, point.Y, 0.0)
        for point in source_points
    ]
    for branch in branches:
        denominator = float(len(branch) - 1)
        for position, node_index in enumerate(branch):
            fraction = float(position) / denominator
            result[node_index].Z = -target_rise * (1.0 - fraction)
    return result


def _solve(points, indexed_edges, anchors, loads, spring_strength):
    goals = List[IGoal]()
    for a, b in indexed_edges:
        line = rg.Line(points[a], points[b])
        goals.Add(
            Spring(
                line,
                line.Length * DEFAULT_REST_LENGTH_FACTOR,
                spring_strength,
            )
        )

    for point, load in zip(points, loads):
        goals.Add(Unary(point, load, 1.0))
    for anchor in anchors:
        goals.Add(Anchor(points[anchor], 1.0e9))

    system = PhysicalSystem()
    for goal in goals:
        system.AssignPIndex(goal, 1.0e-8)
    particle_indices = [
        system.FindParticleIndex(point, 1.0e-8, False)
        for point in points
    ]
    if any(index < 0 for index in particle_indices):
        raise RuntimeError("Kangaroo could not map the network particles.")

    residual = float("inf")
    iterations = 0
    for iterations in range(1, DEFAULT_MAX_ITERATIONS + 1):
        system.Step(goals, True, 1.0e-12)
        residual = system.GetvSum()
        if residual <= DEFAULT_SOLVER_TOLERANCE:
            break
    positions = [
        system.GetPosition(index) for index in particle_indices
    ]
    for anchor in anchors:
        positions[anchor] = points[anchor]
    return positions, iterations, residual


def _coupled_hanging_form(points, indexed_edges, crown, supports,
                          branches, loads):
    mean_load = (
        sum(abs(load.Z) for load in loads) / float(len(loads))
    )
    spring_strength = max(1.0, mean_load * DEFAULT_STIFFNESS_RATIO)
    anchors = list(supports) + [crown]
    total_iterations = 0
    residual = float("inf")
    hanging = None
    apex_clearance = -float("inf")
    monotone = False
    stiffness_runs = 0
    for stiffness_runs in range(1, DEFAULT_MAX_STIFFNESS_RUNS + 1):
        hanging, iterations, residual = _solve(
            points,
            indexed_edges,
            anchors,
            loads,
            spring_strength,
        )
        total_iterations += iterations
        apex_clearance = min(
            point.Z - hanging[crown].Z
            for index, point in enumerate(hanging)
            if index != crown
        )
        monotone = True
        for branch in branches:
            previous_z = hanging[branch[0]].Z
            for index in branch[1:]:
                current_z = hanging[index].Z
                if current_z < previous_z - 1.0e-8:
                    monotone = False
                previous_z = current_z
        if apex_clearance >= -1.0e-8 and monotone:
            break
        spring_strength *= 2.0

    support_z = (
        sum(hanging[index].Z for index in supports)
        / float(len(supports))
    )
    natural_rise = support_z - hanging[crown].Z
    if natural_rise <= 1.0e-9:
        raise RuntimeError("Kangaroo produced no positive hanging rise.")

    return (
        hanging,
        total_iterations,
        residual,
        stiffness_runs,
        natural_rise,
        apex_clearance,
        spring_strength,
        monotone,
    )


def _invert_and_scale(hanging, source_points, crown, supports,
                      natural_rise):
    target_rise = source_points[crown].Z
    if target_rise <= 0.0:
        raise ValueError(
            "The crown node from Stage 02 must be above Z=0."
        )
    support_z = (
        sum(hanging[index].Z for index in supports)
        / float(len(supports))
    )
    scale = target_rise / natural_rise
    thrust = [
        rg.Point3d(
            point.X,
            point.Y,
            (support_z - point.Z) * scale,
        )
        for point in hanging
    ]
    for support in supports:
        thrust[support] = rg.Point3d(
            source_points[support].X,
            source_points[support].Y,
            0.0,
        )
    thrust[crown] = rg.Point3d(
        source_points[crown].X,
        source_points[crown].Y,
        target_rise,
    )
    return thrust, scale


def _branch_polylines(points, branches):
    return [
        rg.PolylineCurve([points[index] for index in branch])
        for branch in branches
    ]


def create_kangaroo_initial_form(nodes, edges, nodal_loads):
    source_points = [rg.Point3d(point) for point in _as_list(nodes)]
    edge_lines = _as_list(edges)
    loads = [rg.Vector3d(load) for load in _as_list(nodal_loads)]
    if len(source_points) < 4:
        raise ValueError("Nodes must contain a radial network.")
    if len(loads) != len(source_points):
        raise ValueError("NodalLoads must contain one vector per node.")
    if any(
        abs(load.X) > 1.0e-9
        or abs(load.Y) > 1.0e-9
        or load.Z > 1.0e-9
        for load in loads
    ):
        raise ValueError(
            "NodalLoads must be downward vertical self-weight vectors."
        )

    (
        indexed_edges,
        crown,
        supports,
        branches,
    ) = _topology(source_points, edge_lines)
    planar = _initial_hanging_points(
        source_points, crown, branches
    )
    (
        hanging,
        iterations,
        residual,
        stiffness_runs,
        natural_rise,
        apex_clearance,
        spring_strength,
        hanging_monotone,
    ) = _coupled_hanging_form(
        planar,
        indexed_edges,
        crown,
        supports,
        branches,
        loads,
    )
    thrust_points, vertical_scale = _invert_and_scale(
        hanging,
        source_points,
        crown,
        supports,
        natural_rise,
    )
    thrust_edges = [
        rg.Line(thrust_points[a], thrust_points[b])
        for a, b in indexed_edges
    ]
    thrust_lines = _branch_polylines(thrust_points, branches)

    monotone = True
    for branch in branches:
        previous_z = thrust_points[branch[0]].Z
        for index in branch[1:]:
            current_z = thrust_points[index].Z
            if current_z > previous_z + 1.0e-7:
                monotone = False
            previous_z = current_z

    total_weight = -sum(load.Z for load in loads)
    converged = (
        residual <= max(DEFAULT_SOLVER_TOLERANCE * 10.0, 1.0e-9)
        and apex_clearance >= -1.0e-8
        and hanging_monotone
        and monotone
    )
    report = (
        "mode=kangaroo_real_selfweight_initial; nodes={0}; edges={1}; "
        "arches={2}; supports={3}; total_weight={4:.6g}kN; "
        "load_distribution=nonuniform; crown_load_applied_once=True; "
        "crown_constraint=geometric_initial_only; "
        "spring_strength={5:.6g}; iterations={6}; residual={7:.6g}; "
        "stiffness_runs={8}; natural_hanging_rise={9:.6g}; "
        "vertical_scale={10:.6g}; apex_clearance={11:.6g}; "
        "crown_is_highest={12}; monotone={13}; converged={14}; "
        "status=initial_form_only_TNA_TNO_required"
    ).format(
        len(source_points),
        len(indexed_edges),
        len(branches),
        supports,
        total_weight,
        spring_strength,
        iterations,
        residual,
        stiffness_runs,
        natural_rise,
        vertical_scale,
        apex_clearance * vertical_scale,
        apex_clearance >= -1.0e-8,
        monotone,
        converged,
    )
    return thrust_edges, thrust_lines, converged, report


try:
    _component_inputs = (Nodes, Edges, NodalLoads)
except NameError:
    _component_inputs = None

if _component_inputs is not None:
    ThrustEdges, ThrustLines, Converged, Report = (
        create_kangaroo_initial_form(*_component_inputs)
    )
