# 图与数据索引

运行 `python run_all.py` 自动生成 8 PNG、12 CSV、1 NPZ 和 1 JSON。

| 图 | 对应数据/记录 |
|---|---|
| bands_transmission.png | transmission.csv；能带由 summary.json 的参数和 profile() 重建 |
| iv_ndr.png | iv_300K.csv、differential_conductance.csv、bias_refinement.csv |
| local_spectra.png | local_spectra.npz：位置、能量、LDOS、单能量 AL/AR、偏压 |
| conservation.png | bond_spectrum.csv；完整能量轴由 local_spectra.npz 给出，R/T 可重算 |
| sensitivity.png | well_sensitivity.csv、temperature_sweep.csv |
| independent_validation.png | continuum_comparison.csv、mesh_convergence.csv、summary.json |
| sampling_failure.png | sampling_failure.csv、barrier_linewidth.csv；局部谱由参数重算 |
| convergence.png | mesh_convergence.csv、quadrature_convergence.csv |

CSV 首行为单位明确的英文列名，电流是 A/m²；图中 I–V 为 A/cm²。能量 eV、空间 nm。NPZ 中 LDOS 为每自旋的纵向局域谱，单位 1/(eV nm)，AL/AR 为每节点每自旋 1/eV；未乘费米占据，不是电荷密度。summary.json 汇总参数、峰谷、验证与窄峰错误。

谱通量无量纲；左右注入分别为 ±T。PVR 取第一处内部峰与其后的第一处谷。窄峰采样错误比较的是固定能窗的透射面积，不是完整电流。不同输出量的物理范围与误差见 docs/MODEL.md 和 docs/VALIDATION.md。
