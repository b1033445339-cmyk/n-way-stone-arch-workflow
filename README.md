# N-Way Stone Rib Arch Workflow

A Rhino 8 / Grasshopper workflow for the parametric design and feedback-driven revision of irregular radial stone rib arches sharing one crown.

Companion material for **A Parametric Workflow for the Design and Feedback Evaluation of Irregular N-Way Stone Rib Arches**, CAADRIA 2027 Contribution #262, by Zhifan Bai, Tao Zhou and Shixing Liu.

The workflow connects parametric geometry, nonuniform self-weight, Kangaroo initial form finding, a compression-only network, section screening, stone segmentation, and joint contact/sliding evaluation. The designer sets inputs, reviews the returned results and chooses the next revision.

## Download and use permission

You may download, keep a local copy, configure inputs and run the tool for your own work. Redistribution, republication, sublicensing and publication of modified versions require separate permission. This public repository uses a custom limited download-and-use permission, not an open-source license. See [LICENSE.txt](LICENSE.txt) for the complete notice.

允许下载、保存本地副本、设置输入参数并运行工具；再分发、重新发布、转授权或发布修改版本须另行取得许可。本仓库公开可获取，采用自定义下载与使用许可，不标称开源项目。完整说明见 [LICENSE.txt](LICENSE.txt)。

## Contents

- `grasshopper/NWayStoneArch.gh`: saved Grasshopper definition.
- `scripts/`: 12 paper-specific Python components, ending at the sliding-deficit envelope.
- `examples/paper_cases.json`: current manuscript inputs H0–H3 and N3/N4/N6, shared study settings and case roles.
- `examples/N6_reference_results.json`: N6 research-mainline reference results, warnings and verification scope.
- `docs/SCRIPT_INTERFACES.md`: input/output descriptions extracted from the distributed scripts.
- `docs/PROVENANCE.md`: verification basis, artifact identity and model scope.
- `SHA256SUMS.json`: file identities for the prepared material.

## Environment and use

Use Rhino 8 with Grasshopper and Python 3, Kangaroo/KangarooSolver, NumPy, SciPy, COMPAS **2.15.1** and COMPAS TNA **0.7.0**. RhinoCommon and the .NET/Kangaroo interfaces are required; these are not standalone command-line Python scripts. Third-party packages and Rhino are not bundled.

1. Prepare the dependencies in Rhino's Python environment. Set the Rhino document units to metres and absolute tolerance to **0.001 m**.
2. Open the saved definition manually in Grasshopper. Inspect the controls and confirm that the required dependencies are available.
3. Set one complete candidate from `examples/paper_cases.json`. Use `N` as the branch count. Feed the same crown radius and top slope to geometry generation and segmentation. Set the shared material, friction and directional-study values explicitly.
4. Run the stages in order: 01 → 02 → 03 → 04 → 05 → 06 → 07A → 07B → 07C → 07D → 07E → 07F. Stage 07A is the file named `07_BlockWeightAndJointAdjacency.py`. Read the component reports before interpreting downstream values.
5. Change a geometric input, then recompute that candidate's own self-weight, network, segmentation and joint results. Keep each result with its own geometry and settings. Record the maximum deficit, its joint/direction, and the separately reported pressure and warnings.

H1 is an example for proceeding through the full paper workflow. H0 intentionally retains a containment failure and therefore has no downstream joint results in the paper. Warnings remain in the completed reference evaluations. The supplementary cases differ in geometry and height, so they do not isolate branch count as the only variable. The manuscript uses N3, N4 and N6 as supplementary cases.

The saved definition and companion sources retain the published implementation. N6 was verified in the research mainline on 26 October 2026; its parameters and selected results are now documented here. This update does not constitute a numerical rerun of the distributed definition or paper-specific script set. See [the provenance note](docs/PROVENANCE.md) for the distinction.

## 中文使用说明

本工具用于不规则多向石砌肋拱的方案生成、评估和修改比较。先在 Rhino 8 中准备 Grasshopper、Python 3、Kangaroo、NumPy、SciPy、COMPAS 2.15.1 和 COMPAS TNA 0.7.0，设置米制和 0.001 m 容差，再由使用者打开 `grasshopper/NWayStoneArch.gh`。

从案例文件选取一整组参数，并显式设置共用材料与方向扫描参数。几何生成与分块使用相同的冠部半径和顶面坡度。每次修改后重新计算该方案自身的自重、受压网络、分块与接缝结果，按阶段报告读取状态。H0为保留的包含检查失败案例；完成计算的案例仍可能有截面警告、接缝张开和正滑移缺额。当前论文以H0–H3为主案例，N3、N4、N6为补充。N6的冠部圆角为0.20 m，完整输入与研究主线参考结果见examples目录。

各参数的取值范围需结合具体设计条件与工程经验初步判断，再评估所选组合。材料参数及简化荷载用于研究比较，工程应用需进一步评估。N6已于2026年10月26日在研究主线完成核验，本次更新整理并公开其参数与必要结果摘要，没有重新运行计算。公开脚本及Grasshopper文件保持原版本，不将研究主线的N6核验等同于这些公开文件已完成N6复测。
