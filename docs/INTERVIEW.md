# 英文简介与复试追问

这份材料用于完成复现后的训练。不要在尚未复现时将工具生成的参考实现称为全部独立编写，也不要将教学参数结果当作实验发现。

## 一分钟英文简介：完成复现后按实际贡献改写

I reproduced a coherent quantum-transport model of a double-barrier resonant tunneling diode. The reference implementation was prepared with AI assistance. My own work focused on [replace this sentence with the derivation, implementation, debugging, or extension you actually completed].

The model uses contact self-energies and a recursive Green-function solver to calculate transmission. I integrated the transverse Fermi supply to obtain current density and examined the first negative differential resistance region. I checked current conservation, compared the solver with a separate banded solve, and studied convergence toward a continuum scattering reference at zero bias.

A useful numerical lesson was that a narrow resonance can produce a large integration error even when the plotted curve looks reasonable. I therefore checked energy spacing, grid offsets, spatial resolution, and bias resolution separately. The potential is prescribed, so this is not a self-consistent NEGF–Poisson simulation or an experimentally calibrated device model.

## 关键追问

1. **为什么从上一个封闭量子问题转向开放系统？** 两个电子库有不同化学势，需要注入和逃逸，不能用有限盒子本征态直接代替透射问题。开放量子系统仍可有准束缚共振。
2. **自能为什么是复数？** 半无限接触消去后给端点带来能量偏移与传播逃逸；推迟支虚部非正。带外没有传播通道，虚部为零，不可套常数展宽。
3. **为什么不用一个很小的 iη？** 人为吸收会引入额外损失，改变幺正性和窄峰。这里用接触解析推迟极限，没有额外 η。
4. **如何判断接触匹配正确？** 均匀链 T=1、R≈0；自能满足二次方程和因果支；连续极限与散射参考相合。单一守恒检查不足以证明参数正确。
5. **峰顶接近 1 是否意味着电流很大？** 不一定。厚势垒可保留对称结构峰顶，却使峰窄；电流依赖整个透射谱与供给权重的积分。
6. **负微分电阻为什么不违反被动性？** 在约定方向下 J 与 V 同号，直流功率仍为正；dJ/dV<0 表示增量斜率，不等于负的总直流电阻。
7. **为什么更大偏压反而减小电流？** 共振能级的位移、左接触传播阈值与费米供给共同作用；不能只用化学势窗口变大推断电流单调上升。
8. **一维结果如何变成 A/cm²？** 横向连续动量积分提供每面积通道数，自旋只计一次；eV 的供给与能量积分各有一个 q 换算，最终系数 q³m/(2π²ħ³)。再除 1e4 转为 A/cm²。
9. **LDOS 是载流子浓度吗？** 不是。这里是每自旋纵向局域谱，需乘接触费米占据并做能量和横向积分才得到三维密度。本项目未计算自洽电荷。
10. **为什么起点改变会影响谱面积？** 峰宽小于采样步长时，不同网格可能恰好命中峰顶或完全错过。配对求积也可能共同漏峰，需要分辨率和独立参考。
11. **如何区分三种收敛？** 空间网格改变 Hamiltonian；能量积分改变谱加权精度；偏压网格改变极值定位。三者互不替代，且都不是材料模型验证。
12. **连续参考为何只比较零偏压？** 当前级联参考使用严格分段常数势；有限偏压为分段线性背景，需要 Airy 或单独控制分层逼近。不能把零偏压一致冒充全部偏压独立验证。
13. **PVR 如何提取？** 第一处内部局部最大、随后第一处局部最小；另外细化局部偏压。不能把高偏压背景上升的末端当作第一共振峰。
14. **温度变化能说明什么？** 当前只说明电子库热展宽；没有声子、自热和温度依赖材料参数，不能声称解释了全部实验温度效应。
15. **厚度敏感性是否就是工艺良率？** 不是。这里是三种确定性尺寸情景；真实良率需尺寸分布、相关性、界面粗糙与工艺标定。改变阱宽还改变给定电压下降长度，要讨论控制变量。
16. **为什么没有把上一个项目的 Poisson 直接拼进来？** 非平衡占据需要左右接触谱、正确横向积分与电荷参考，电中性和边界条件也需验证；直接拼接可能产生假自洽和假滞回。
17. **下一步最有信息量的扩展是什么？** 先引入左右不对称势垒，预测共振峰顶下降并检验互易性；随后才增加变质量和电荷反馈。每步保留独立证据。
18. **你本人做了什么？** 按事实区分参考代码、独立推导、自己重写的模块、发现的错误、验证结果与新扩展。若仅运行过脚本，就说复现；不要说独立开发了完整 TCAD 平台。

## 展示时建议准备的三张图

先用 I–V 引出问题，再用偏压透射谱解释对齐，最后用网格起点导致的谱面积错误说明验证意识。能按问题—证据—局限讲清楚，比背下所有小数更有效。
