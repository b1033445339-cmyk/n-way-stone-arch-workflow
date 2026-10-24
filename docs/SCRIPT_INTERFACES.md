# Script interfaces

The following interface notes are extracted verbatim from the distributed source files. Geometry is in metres, forces in kN and pressures in MPa. Item/List access is significant when constructing Python 3 components manually.

## 01_ParametricGeometry.py

```text
Stage 01 - deterministic parametric geometry for N radial half arches.

Grasshopper Python inputs:
    N                  int          Number of half arches / supports (>= 3)
    H                  float        Overall crown top height above Z=0
    RotationAngles     list[float]  One plan angle in degrees per half arch
    Spans              list[float]  One crown-to-support span per half arch
    FootWidths         list[float]  One transverse width per half arch
    CrownDepth         float        Minimum vertical depth at the crown
    CrownFilletRadius  float        Plan fillet radius at crown intersections
    TopSlopePercent    float        Optional outward top fall, 0 to 5 percent
    CrownBlockRadius   float        Shared crown-zone / first-joint control

Use Item Access for N, H, CrownDepth, CrownFilletRadius, TopSlopePercent and
CrownBlockRadius. Use List Access for RotationAngles, Spans and FootWidths.
Empty RotationAngles gives an even 360-degree distribution. One Span or
FootWidth is broadcast to all arches. Unconnected crown inputs use their
DEFAULT values. An unconnected TopSlopePercent preserves the historical
horizontal top (0 percent); an unconnected CrownBlockRadius uses the same
1.60 m default as Stage 06.

Outputs:
    WholeSolid, CenterLines, SupportPoints, Valid, Report

This component performs geometry generation only. It does not use Kangaroo and
does not claim structural equilibrium. Every half arch is a deterministic
radial wall with a quadratic lower arch edge, an exact support location and a
constant transverse FootWidth. At zero slope all branches retain the historical
common horizontal top plane at Z=H. The wall is deepest at the support and has
the internal CrownDepth at the common crown. CrownFilletRadius rounds the
actual plan intersections of the half arches; it never adds a sphere or
ellipsoid. Block segmentation belongs to a later stage. A positive
TopSlopePercent retains a horizontal central cap of CrownBlockRadius and then
lowers each branch top radially outward at the requested gradient; it does not
cut the later crown block. The same CrownBlockRadius value is intended to feed
Stage 06 in parallel; no Stage 06 output is fed back into Stage 01.
```

## 02_RealSelfWeight.py

```text
Stage 02 - lump the real solid self-weight onto a radial thrust network.

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
```

## 03_KangarooInitialForm.py

```text
Stage 03 - Kangaroo initial thrust form under real nonuniform self-weight.

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
```

## 04_TNAEnvelopeCheck.py

```text
Stage 04 - compression-only TNA equilibrium and solid containment check.

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
```

## 05_SectionKernCheck.py

```text
Stage 05 - no-tension section-kern check for the TNA network.

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
```

## 06_UtilisationGuidedSegmentation.py

```text
Stage 06 - utilisation-guided radially graded stone segmentation.

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
```

## 07B_JointResultants.py

```text
Stage 07B - unified ideal and horizontal-study joint resultants.

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
```

## 07C_JointOpeningAndPressure.py

```text
Stage 07C - actual-joint eccentricity, opening and contact pressure.

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
```

## 07D_PartialContactPressure.py

```text
Stage 07D - compression-only partial-contact pressure redistribution.

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
```

## 07E_SlidingIndicator.py

```text
Stage 07E - directional sliding indicator, without contact-strength verdict.

Inputs: JointPlanes, JointIDs, CaseAngles, CaseJointIDs, ContactResultants,
MaxPressures (MPa, validity check/report only), FrictionCoefficient,
SlidingSafetyFactor. Lists are case-major in the exact JointIDs order.
Outputs: SlidingUtilisations, SlidingCaseKeys, Valid, Report.
eta_s = SlidingSafetyFactor*T/(FrictionCoefficient*N).
N is the positive projection magnitude supplied by Stage 07B, not a signed
tension check. Pressure is reported separately; no measured strength,
crushing assessment, connector quantity or structural safety is inferred.
```

## 07F_SlidingDeficit.py

```text
Stage 07F - per-joint sliding-deficit envelope (paper version).

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
```

## 07_BlockWeightAndJointAdjacency.py

```text
Stage 07A - exact block self-weight and actual joint adjacency.

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
```
