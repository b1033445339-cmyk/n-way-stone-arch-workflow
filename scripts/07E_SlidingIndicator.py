#! python 3
"""Stage 07E - directional sliding indicator, without contact-strength verdict.

Inputs: JointPlanes, JointIDs, CaseAngles, CaseJointIDs, ContactResultants,
MaxPressures (MPa, validity check/report only), FrictionCoefficient,
SlidingSafetyFactor. Lists are case-major in the exact JointIDs order.
Outputs: SlidingUtilisations, SlidingCaseKeys, Valid, Report.
eta_s = SlidingSafetyFactor*T/(FrictionCoefficient*N).
N is the positive projection magnitude supplied by Stage 07B, not a signed
tension check. Pressure is reported separately; no measured strength,
crushing assessment, connector quantity or structural safety is inferred.
"""
import math
import Rhino.Geometry as rg
UTILISATION_TOLERANCE=1.0e-9
NORMAL_FORCE_EPSILON=1.0e-8

def _as_list(value):
    if value is None:return []
    if isinstance(value,str):return [value]
    try:return list(value)
    except TypeError:return [value]

def _unwrap(value):
    return value.Value if hasattr(value,'Value') else value

def _positive(value,name):
    number=float(value)
    if not math.isfinite(number) or number<=0:raise ValueError(name+' must be finite and positive.')
    return number

def _case_key(index,angle,identifier):
    return 'case[{0}]/angle[{1:.6g}]/{2}'.format(index,angle,identifier)

def create_directional_sliding_indicator(joint_planes,joint_ids,case_angles,case_joint_ids,contact_resultants,max_pressures,friction_coefficient,sliding_safety_factor):
    planes=[_unwrap(x) for x in _as_list(joint_planes)]
    ids=[str(x) for x in _as_list(joint_ids)]
    angles=[float(x) for x in _as_list(case_angles)]
    row_ids=[str(x) for x in _as_list(case_joint_ids)]
    forces=[_unwrap(x) for x in _as_list(contact_resultants)]
    pressures=[float(x) for x in _as_list(max_pressures)]
    mu=_positive(friction_coefficient,'FrictionCoefficient')
    gamma=_positive(sliding_safety_factor,'SlidingSafetyFactor')
    n=len(ids); expected=n*len(angles)
    if not n or not angles or len(set(ids))!=n or len(planes)!=n or any(len(a)!=expected for a in [row_ids,forces,pressures]):
        raise ValueError('Directional inputs require a complete JointIDs-sized block per direction.')
    if any(not math.isfinite(a) for a in angles):raise ValueError('CaseAngles must be finite.')
    normals=[]
    for p in planes:
        if not isinstance(p,rg.Plane):raise ValueError('Invalid JointPlane.')
        normal=rg.Vector3d(p.Normal)
        if not normal.Unitize():raise ValueError('Invalid joint normal.')
        normals.append(normal)
    values=[];keys=[];failures=[]
    for case,angle in enumerate(angles):
        start=case*n
        if row_ids[start:start+n]!=ids:
            failures.append((case,angle,'case_joint_id_order_mismatch'));continue
        for j,identifier in enumerate(ids):
            i=start+j
            try:
                force=forces[i];p=pressures[i]
                if not isinstance(force,rg.Vector3d):raise ValueError('Invalid contact resultant.')
                if not math.isfinite(p) or p < -1e-12:raise ValueError('Peak pressure must be finite and non-negative.')
                N=rg.Vector3d.Multiply(force,normals[j])
                if not math.isfinite(N) or N<=NORMAL_FORCE_EPSILON:raise ValueError('Normal magnitude is below the admissible gate.')
                T=(force-normals[j]*N).Length
                if not math.isfinite(T):raise ValueError('Invalid tangential magnitude.')
                # Preserve the production operation order and screening tolerance.
                required_friction=gamma*T/N
                eta=required_friction/mu
                values.append(eta)
                if eta>1.0+UTILISATION_TOLERANCE:keys.append(_case_key(case,angle,identifier))
            except Exception as e:failures.append((case,identifier,str(e)))
    valid=not failures and len(values)==expected
    if not valid:
        # Never publish a silently filtered subset of directional results.
        values=[];keys=[]
    report=('stage=directional_sliding_indicator; angle_cases={0}; joints_per_case={1}; '
            'flattened_rows={2}; sliding_case_rows={3}; friction_coefficient={4}; '
            'prescribed_sliding_factor={5}; pressure_max_MPa={6}; failures={7}; valid={8}; '
            'N_is_projection_magnitude_not_signed_audit; pressure_assessed_separately; '
            'no_contact_strength_or_connector_capacity_verdict; research_scenario_only').format(
                len(angles),n,len(values),len(keys),mu,gamma,max(pressures) if pressures else None,failures,valid)
    return values,keys,valid,report

try:
    _component_inputs=(JointPlanes,JointIDs,CaseAngles,CaseJointIDs,ContactResultants,MaxPressures,FrictionCoefficient,SlidingSafetyFactor)
except NameError:_component_inputs=None
if _component_inputs is not None:
    SlidingUtilisations,SlidingCaseKeys,Valid,Report=create_directional_sliding_indicator(*_component_inputs)
