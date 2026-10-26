# Verification basis and model scope

## Current contents

The repository contains 12 paper-specific component scripts, one saved Grasshopper definition, the current manuscript's H0–H3 and N3/N4/N6 parameters, and the N6 reference-results summary. H0 is an intentional sampled-containment failure in the main height-adjustment sequence; it has no downstream joint results. File identities are listed in `SHA256SUMS.json`.

## Numerical evidence

The retained verification records for the paper workflow cover analytical contact examples, an ideal zero-deficit case, an invalid directional case and the main/supplementary cases available before the N6 addition. The final H3 record reports 32 blocks, 31 joints, 240 directions, a maximum sliding deficit of 14.596892268097193 kN and maximum joint pressure of 0.6523723590797497 MPa.

N6 was verified on 26 October 2026 using the research-mainline candidate `attempt10`. Its full Stage 01–07H execution reported `COMPLETED` and `EXECUTION_PASS=True`, followed by independent saved-data checks and a reopened geometry artifact. The source evidence hashes were checked when preparing this update; no new numerical analysis was run. `examples/N6_reference_results.json` contains the selected results, warnings and verification scope through the paper's joint-assessment stages.

N6 has 36 blocks and 35 joints, with 240 directions and 8,400 joint-direction rows. Its maximum sliding deficit is 6.4120687296 kN at arch[3]/joint[5], 228 degrees; its maximum joint pressure is 1.4563088896 MPa at arch[6]/joint[6], 290 degrees. Its 0.20 m crown fillet differs from the H-series value of 0.30 m. Section warnings, opening and positive sliding deficit remain. The first fully completed candidate was retained, rather than selecting a candidate by pressure or deficit magnitude.

## Version and reproduction scope

The research-mainline source set used for N6, the distributed paper-specific scripts and the saved Grasshopper definition are separately versioned artifacts. The distributed Stage 04, 05 and 06 sources differ from the current research-mainline files; distributed 07E/07F are the paper-specific sliding-indicator/deficit components. This publication update does not constitute an N6 numerical rerun using the distributed definition or every distributed component wrapper. The scripts, definition and interfaces are unchanged by this data/documentation update.

When manually rebuilding components, follow `docs/SCRIPT_INTERFACES.md` and use one candidate's complete geometry and settings throughout the pipeline. The published components end at the sliding-deficit envelope; empirical dowel quantities and display geometry are outside this repository's scope.

## Interpretation

The numerical model selects a compression network under stated constraints and evaluates sampled containment, sections and actual joints. Component completion and containment are separate from structural safety. The horizontal study is a pseudo-static directional comparison. Material and friction inputs are assumed research values. Parameter feasibility, block equilibrium, material response and dynamic behaviour require further investigation, as discussed in the paper.

## 中文说明

本仓库仅提供当前论文采用的H0–H3、N3/N4/N6参数。H0是主案例序列中保留的包含检查失败方案，并非备用历史案例。N6参数与结果来自已核验的研究主线；本次仅整理发布并核对既有证据。公开论文脚本、Grasshopper定义与研究主线分别版本化，不能把主线核验直接表述为公开文件已完成N6复测。公开组件止于滑移缺额包络，计算完成不等于结构安全。
