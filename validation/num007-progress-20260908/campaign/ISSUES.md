# 首轮问题台账（尚未集中修复）

> 当前裁决：以末尾“源码设计复审与台账裁决检查点”为准；本文保留历史观察，旧标题和旧严重程度不再单独作为最新判错依据。

构建：PR #3 `bdc3ece61b3847a782873bfd2654521c9808efb5`。当前没有完整科学通过结论。原库 273 项分别记录在 `existing_cases.jsonl` / `reference_cases.jsonl`；同一路径以最后一条有效记录为准。

## 已复现的软件差异

### EXP-001 — 实际 GUI 导出绕过当前化学/电子态数据构建

- 案例：OLD-001，`gaussian_validation_set/01_H3plus_3c2e.fch`；OLD-002，`gaussian_validation_set/02_XeF2_3c4e.fch`。
- 严重度：Major；类别：C，显示/导出及 provenance。
- 操作：当前 COV 载入文件，紧凑多中心 MO 图，分别选择两个 E′ 成员，点击 Export diagram + metadata。
- 屏幕：紧凑图含两个能级行、三个原始 MO；MO1 为 A1′、占据 2，MO2/3 为 E′、各占据 0。电子态 +1 / singlet，producer FCHK 和同源日志可用。
- 实际导出：`mode=valence-central`，7 个显示轨道、7 行；`wavefunction_source=Unknown`，charge/multiplicity/alpha/beta 为 null。不符合屏幕和输入。
- 证据：`evidence/OLD-001-structure.png`、`OLD-001-homo.png`、`OLD-001-mo1-tooltip.png`、`OLD-001-mo2a-tooltip.png`、`OLD-001-mo2b-tooltip.png`、`OLD-001-ocr.json`；实际导出包 `OLD-001-original.mo.{png,svg,json,csv}`。
- 代码依据：CMake 对 `mo_diagram_v2.cpp` 整个编译单元把 `build_mo_diagram_data` 宏改名为 legacy；part8 的 export 函数内部同名调用因此进入 legacy；屏幕进入 `mo_diagram_chemistry.cpp` 的当前公开实现，后者补电子态并执行化学筛选。
- 正确结果：同一输入/设置下，屏幕与四种导出的 raw MO 集合、行、占据、模式、能量轴和 provenance 一致，不以 null/Unknown 覆盖已知 producer 数据。
- OLD-002 复现：紧凑屏幕为原始 MO 6/11/14 三行；实际导出 21 行、mode=valence-central，电子态仍 Unknown/null。证据为 `OLD-002-mo11-tooltip.png`、`OLD-002-mo6-tooltip.png`、`OLD-002-ocr.json` 及 `OLD-002-original.mo.{png,svg,json,csv}`。两条 Xe–F 实线连接在 `OLD-002-structure.png` 的 opacity=0.02 下已确认。
- 待集中修复：使导出明确调用公开最终构建入口或复用同一已构建数据。增加真实 H3+、UHF、π-only、配体场和没有化学数据的负例；导出测试必须比较公开 builder 而非 legacy 自洽。

### NUM-001 — 球谐 AO 的奇数 |m| 相位约定与 Gaussian 不匹配

- 严重度：Critical；类别：A，AO/MO 实空间数学及可能的派生分析。
- 独立来源：Gaussian 16 `cubegen.exe` 读取同一 FCHK，AMO/BMO 指定相同自旋轨道，Bohr 轴向网格；通过真实 `CudaOrbitalEvaluator` 和 OpenGL texture readback 取 COV 网格。不是根据截图猜测。
- 初始样本：OLD-003/004/005/006，最大 sampled NRMS 分别约 1.048、0.784、0.122、1.96；同文件也有误差约 1e-6 的轨道。整体正负相位对齐后仍有差异。
- 定位：`orbital.cu::real_solid_harmonic` 的纯 d/f/g 奇数 |m| 多项式带负号，而 FCHK pure shell 系数以 identity/sign=+1 导入。CPU `orbital_symmetry.cpp` 也有相同多项式约定，需审查所有消费者，不能只改渲染。
- 诊断实验（不算正式修复通过）：测试程序仅在内存副本中翻转纯壳奇数 |m| 系数，四例全部 sampled MO 达到 NRMS <= 6.86e-6。结果在 `phase_diagnostic.json`；正式程序和默认测试仍未修复。
- 原始失败 cube：`evidence/OLD-xxx-grid-failure/`；全库扫描记录 `grid_cases.jsonl`。通过案例的 disposable cube 用后移除，仅保留可复算的网格/命令/数值；失败参考保留。
- 科学依据：Gaussian cubegen 官方手册转载原文 https://theochem.mercer.edu/chm295/g09ur/u_cubegen.htm （Gaussian 主站 502）；IOData FCHK source https://iodata.readthedocs.io/en/stable/_modules/iodata/formats/fchk.html 及 HORTON basis conventions https://theochem.github.io/horton/2.1.1/tech_ref_gaussian_basis.html （2026-09-05 查阅）。Gaussian pure conventions 与 HORTON 一致，后者不采用 Condon–Shortley 相位。
- 待集中修复：明确统一内部 AO 约定；密度双指标、overlap、CPU symmetry/chemistry、CUDA、Molden 兼容一起核对。增加纯 d/f/g、Cartesian/SP、混合壳、α/β、非轴对齐分子，以及错误相位负例。不能用全局红蓝互换解释局部 AO 相位错配。
- 当前范围：全库 273 文件扫描完成，生产基线 26 通过、247 不通过。对全部 247 个失败文件的原始失败 cube 进行内存相位诊断，247 个均达到上述容差，记录在 `phase_diagnostic_cases.jsonl`。这不算生产修复通过；这里只验证被测网格/轨道，不能将 representative sampling 宣称每一个轨道或空间积分均已通过。

### CHEM-001 — 成键分布汇总携带默认 UND 权重，百分比超过 100%

- 严重度：Moderate；类别：B/C，派生分布及显示。
- 案例：OLD-003，选择 source MO 26 / internal 25，T2g，2e，-0.0344978169 Ha。
- 截图：`evidence/OLD-003-classification-percent.png`；OCR 在 `OLD-003-ocr.json`，已回看画面校正。界面为 bonding 78%、antibonding 0%、nonbonding 22%、UND 503%。
- 原因：`include/cov/model.hpp:106` 的 `OrbitalBondingDistribution` 默认 undetermined=1；`src/analysis/orbital_chemistry.cpp:862` 的 aggregate_bonding 未将累加器清零，就加权累加各 interaction，最后所有分量除 total。于是未知分量被额外增加 1/total。
- 正确语义：有有效 interaction 时，从零开始汇总四个互斥分布分量并归一化；无有效 interaction 时单独返回 UND=1 的降级状态。不应只截断 UI 百分比或随意改全部默认对象语义。
- 待集中修复/回归：空列表、全 NotApplicable、单项、混合成键/反键/非键/未知、缩放所有 interaction 权重均须检查非负、有限、总和为 1 和缩放不变；真实 FCHK 全 MO 扫描分布不变量，并比较 UI 与数据。

### OLD-003 的显示采集补充

- 六条 Ti–F coordination contacts 在 opacity=0.02 下齐全；正常 MO 使用 opacity=0.92 / isovalue=0.0300 / grid=128。
- 紧凑 8 行；切换展开 17 行；恢复紧凑仍为 8 行。T2g* source MO 42 与 43 可独立选取，T2g source MO 26 可选，不会因选择额外删除行。
- 分裂符号 Δπ 置于侧边；tooltip 为 π-donor，0.369541 Ha、pair confidence=0.909，与 0.335043428-(-0.0344978169) 的数值一致。非线性轴不代表线性比例。
- 导出再次复现 EXP-001：27 行、valence-central、Unknown/null 电子态，屏幕为 ligand-field 8 行；实际包 `OLD-003-original.mo.{png,svg,json,csv}` 已归档，没有覆盖旧输入。
- 结构、完整相关图、逐成员 tooltip、切换前后分别为 `OLD-003-structure.png`、`OLD-003-compact-complete.png`、`OLD-003-t2g-star-b.png`、`OLD-003-t2g-star-c.png`、`OLD-003-t2g-b.png`、`OLD-003-expanded-upper.png`、`OLD-003-expanded-complete.png`、`OLD-003-compact-restored.png`；载入身份 `OLD-003-loaded.png`；轴模式 `OLD-003-compact-upper.png`。
- 判定：不通过（NUM-001、EXP-001、CHEM-001）；其余配体场化学筛选和 REG-001 断言仍待独立核验，不把交互稳定当作科学通过。

### OLD-004 的显示采集补充

- `04_ZnCl4_2minus_Td.fch`：5 原子、103 basis/MO、-2 / singlet、50/50 电子；同源日志附加 TD 对称性。最低 opacity 的四条 Zn–Cl 联系齐全，需重置继承的相机后观察，不能把视角遮挡当缺键。
- 紧凑 5 行、展开 15 行、恢复紧凑 5 行；T2 source MO 46/47 和 E source MO 37 均可独立选择，选择前后未复现额外消失。E 行为 -0.160095621 Ha，近似非键，Δπ=0.003720 Ha，pair confidence=0.814。
- T2 source 47 为 +0.0490927801 Ha / 2e，E source 37 的 Zn 3d 权重 96.6%。正占据轨道能量是输入事实，不单独作为 COV 出错的依据；该气相二价阴离子的电子束缚、弥散基组和稳定性需另行质控。
- 同一个 E MO 的 tooltip 显示 nonbonding / 81%，下方 Selected MO chemistry 却为 mixed (UND)，需在集中修复中核对原始 chemistry 与配体场派生标注的来源和一致性。
- 实际导出再次复现 EXP-001：27 行、valence-central、Unknown/null 电子态，而屏幕紧凑 5 行。实际包为 `OLD-004-original.mo.{png,svg,json,csv}`；仅移动本轮刚导出的四个文件到证据目录，没有覆盖输入或旧证据。
- 证据：`OLD-004-loaded.png`、`OLD-004-structure.png`、`OLD-004-mo-normal.png`、`OLD-004-compact-complete.png`、`OLD-004-t2-b.png`、`OLD-004-t2-c.png`、`OLD-004-e-a.png`、展开完整图及恢复图、`OLD-004-ocr.json`。判定不通过（NUM-001、EXP-001），其他化学分析与输入质控仍分项保留。

### OLD-005 的显示采集补充

- 苯 `05_benzene_D6h_def2SVP.fch`：12 原子、114 basis/MO、0 / singlet、21/21；producer D6H 标签可用。当前文件为固定几何单点，不宣称是已验证极小。
- 最低 opacity=0.02 下六条 C–H 实线和环边连接可见，没有穿环中心的跨环伪键；相机角度影响遮挡，正常 MO 证据使用 opacity=0.92 / isovalue=0.03 / grid=128。
- 紧凑 π-only 保留 source 17、20、21、22、23、29，6 个 MO / 4 行 / 6 个电子；20/21 同为 -0.266338021 Ha、E1g、占据各 2，可以独立选择；29 为 +0.158962250 Ha、B2g、占据 0。展开 19 行、恢复紧凑 4 行，选中 29 不丢失。
- CHEM-001 复现：MO20 Selected chemistry 为 bonding 64%、antibonding 32%、nonbonding 4%、UND 175%；MO29 为 21%、79%、0%、UND 21%。不能将 tooltip 的子空间分类与原始 interaction 汇总百分比混作同一定义。
- EXP-001 复现：当前 π-only 屏幕 4 行，实际导出却为 25 行、valence-central、Unknown/null 电子态，且丢失部分 π 成员与当前选择。实际包 `OLD-005-original.mo.{png,svg,json,csv}`；只移动本轮新导出到证据目录。
- 证据：`OLD-005-loaded.png`、`OLD-005-structure.png`、`OLD-005-mo-normal.png`、`OLD-005-compact-complete.png`、`OLD-005-homo-a.png`、`OLD-005-homo-b.png`、`OLD-005-pi-antibonding.png`、展开完整图/恢复图/导出状态、`OLD-005-ocr.json`。OCR 已回看，数值以 FCHK 和可读截图校正。判定不通过（NUM-001、EXP-001、CHEM-001）。

### NUM-002 — 截断 MO 块导致 overlap/Mayer 不可用

- OLD-020 `cp_anion_pi_vmd.fch` 为 400 AO / 396 个保留 MO，不是方阵。COV 返回空 AO overlap 和空 Mayer 列表。
- 独立 IOData/GBasis 直接对原始基组积分，S 最小特征值 2.02227e-7，条件数 1.39216e8；CᵀSC 最大非对角偏差 1.02922e-6，Tr(PS)=35.9999999944。因此不能把输入直接判为无效或非键。
- 由不完整 C 唯一反推完整 S 在数学上做不到；原始基组信息足够解析积分。待集中处理：增加受控的解析 AO overlap 路径，明确病态基组容差与 provenance；未得到 S 时必须如实 unavailable，不能给确定化学结论。具体降级显示仍待该例 UI 核对。
- 证据：`overlap_cases.jsonl`、`analysis_comparison_cases.jsonl` schema_version=2；Mayer 独立环边约 1.1377–1.1378、C–H 约 1.2081；不把小的负远程 Mayer 当成应画的化学键。

### SYM-001 — 偶数阶 Dnh 使用错误的撇号标签并漏掉 B 类

- 严重度：Major；类别：B，对称性派生；OLD-006 旋转且 NoSymm 的苯。几何识别 D6h，但 source MO20/21 的 HOMO 被标 E1''，source MO29 的最高 π 反键轨道为 N/A。
- 正确 D6h π 子空间标签为 A2u、E1g、E2u、B2g；同几何刚体旋转不能将 g/u 变为撇号。MO20/21 分别 -0.266338011/-0.266337921 Ha、各占据 2，29 为 +0.158961980 Ha、占据 0。不能把这类缺标签归因于 FCHK 不直接提供标签就停止分析。
- 源码 `src/analysis/orbital_symmetry.cpp::classify_dnh` 无条件用水平镜面生成 ' / '' 后缀，且一维 Cn character<0 时直接 unavailable；因此偶数 n 的 B1/B2 缺失。集中修复需要偶/奇 n 分开，以反演 character 确定偶数 n 的 g/u，并确定参考 C2 轴约定；不能替换苯的固定 MO 编号。
- 独立参考（2026-09-05 查阅）：Purdue University benzene MO https://www.chem.purdue.edu/gchelp/orbs/benzene.html ；Newcastle University / Jon Goss D6h character table https://www.staff.ncl.ac.uk/j.p.goss/symmetry/D6h.html 。前者为定性标签参考，精确数值只对照本地输入。
- UI 证据 `OLD-006-loaded.png`、`OLD-006-homo-a.png`、`OLD-006-homo-b.png`、`OLD-006-pi-antibonding.png` 和 `OLD-006-ocr.json`；原 OCR 字符有误读，以已查看的实际图和数据校正。
- 结构最低 opacity=0.02 无跨环伪键；紧凑仍为 4 行/6 MO，展开 19 行，恢复 4 行且选择不丢失。实际导出再次复现 EXP-001：25 行、valence-central、Unknown/null；四文件 `OLD-006-original.mo.*` 为本轮实际 GUI 导出，未改旧输入。
- 回归应包括 D3h/D5h 不变、D4h/D6h 各一维与二维表示、随机旋转、原子重排、无 producer 标签和低对称性负例；并与 NUM-001 的 CPU AO 相位修复联合复验。

### OLD-007 较大基组对照

- 苯 `07_benzene_D6h_def2TZVP.fch`：12 原子、90 壳、222 AO/MO、0/singlet、21/21，producer D6H 标签可用。紧凑 4 行/6 MO（17、20、21、22、23、30）/6e，展开实测 16 行，恢复 4 行；不能沿用 def2SVP 的 19 行记录。
- MO21 独立选择成功，E1g、-0.266970103 Ha、占据 2；MO30 为 B2g、+0.149339457 Ha、占据 0，π 子空间 92%。与 NoSymm 苯比较支持 SYM-001 位于 derived 分支；正常标签显示通路未统一损坏。
- 结构 opacity 0.02 下环边及 C–H 可见，无跨环中心伪键。正常 MO 用 opacity 0.92 / isovalue 0.03 / grid 128；更换成员不改变紧凑选择集合。正常表面仍受 NUM-001 的纯球谐相位错误影响，不计作科学通过。
- EXP-001 再现：实际导出 26 MO/26 行、valence-central、Unknown/null 电子态，与屏幕 4 行不一致。四文件 `OLD-007-original.mo.*` 是本轮 GUI 新导出后移入证据目录，未覆盖原输入。
- 证据 `OLD-007-loaded.png`、`OLD-007-structure.png`、`OLD-007-mo-normal.png`、`OLD-007-compact-complete.png`、`OLD-007-homo-b.png`、`OLD-007-pi-antibonding.png`、`OLD-007-expanded.png`、`OLD-007-compact-restored.png`、`OLD-007-export-state.png` 和 `OLD-007-ocr.json`。OCR 错读 21/20、撇号等，以实际截图和原数据校正。

### SYM-002 — 低对称点群缺少可计算的 MO 标签推导

- 严重度：Major；类别：B。OLD-008 反式丁二烯 `08_trans_butadiene_C2h_def2SVP.fch`，几何 C2h、NoSymm 单点。86 个 MO 的标签全为 unavailable；紧凑图四个 π MO 同样 N/A。
- 独立检验：`audit_symmetry.py` 用 IOData/GBasis 直接读取原始基组和 MO 系数，在 768 个确定性采样点计算反演、分子平面反射及法向 C2 操作的表示；先做同元素几何置换检查，再检查同自旋能量子空间的最小二乘残差。不使用 COV 的 AO evaluator 或标签推导代码作为真值。
- OLD-008 四个 π 轨道 source 14/15/16/17 的能量分别为 -0.338438876、-0.243839600、-0.0267780274、+0.0703011038 Ha，表示为 Au/Bg/Au/Bg，最大相对采样残差 4.1917e-14；全部 86 MO 可赋标签。OLD-038 的同类 C2h 输入也完成 252 MO 独立标签检验，UI 尚待进行。逐表示矩阵、残差和几何容差保存在 `symmetry_reference_cases.jsonl`。
- C2h 表示解释使用已识别的点群；独立步骤验证具体几何操作及波函数变换，不宣称另行独立完成所有点群识别。
- 源码 `src/analysis/orbital_symmetry.cpp` 的 derive 分支仅覆盖 Td、Oh 和 n>=3 的 Dnh；其他点群在缺 producer 标签时无推导路径。应使用统一的对称操作/角色投影和残差门槛补齐，不能把 N/A 按固定分子/MO 编号硬替换。
- 独立方法参考：GBasis 官方 MO evaluation https://gbasis.qcdevs.org/tutorial/Evaluations_basis_and_potential.html ；Newcastle University/Jon Goss C2h character table https://www.staff.ncl.ac.uk/j.p.goss/symmetry/C2h.html （2026-09-05）。
- UI：10 原子、42 壳、86 basis/MO、0/singlet、15/15；低 opacity=0.02 见三条 C–C 和六条 C–H 实线，无远程伪键；正常 MO 为 0.92/0.03/grid128。紧凑四个 π MO/4e、展开 21 行、恢复四行且 source16 选择不丢失。HOMO15 与 LUMO16 均独立选择，tooltip 95% π 分类；其空间表面仍受 NUM-001 影响。
- 实际 GUI 导出再次复现 EXP-001：21 行、valence-central，而屏幕为四行 π-only。四文件 `OLD-008-original.mo.*` 已从本轮新导出移动归档，不覆盖原始计算。
- 证据：`OLD-008-loaded.png`、`OLD-008-structure.png`、`OLD-008-mo-normal.png`、`OLD-008-homo.png`、`OLD-008-lumo.png`、`OLD-008-compact-complete.png`、`OLD-008-expanded.png`、`OLD-008-compact-restored.png`、`OLD-008-export-state.png`、`OLD-008-ocr.json`。OCR 已回看，字符错读不作为数值依据。判定不通过；原始单点不宣称已验证极小。

### OLD-009 乙烯及 UI-001 — 非配合物误用配体场摘要措辞

- OLD-009：6 原子、24 壳、48 AO/MO、0/singlet、8/8，producer D2H；最低 opacity 0.02 检查四条 C–H 和一条 C–C 均为正常实线，无错误虚线。正常 MO 为 opacity 0.92 / isovalue 0.03 / grid 128。
- source8 HOMO 为 B3u、-0.285158906 Ha、2e；source9 LUMO 为 B2g、+0.012050992 Ha、0e，均可独立选择。由原基组和 MO 系数在输入 xyz 坐标系下独立计算 D2h 操作，两个标签一致，最大残差 1.052e-15。一个二中心 π 通道、两个 MO、2e；普通孤立 C=C 不自动要求多环/多共轭的 π-only 视图。
- 当前 valence-central 图为 11 行，展开/恢复仍为 11 行，选择不丢失；没有可额外隐藏的框架不属于缺陷。实际 GUI 导出却为 16 行/16 MO，电子态 Unknown/null，继续复现 EXP-001；四文件已归档为 `OLD-009-original.mo.*`。
- UI-001（Minor）：非金属乙烯的普通价层摘要写成 `ligand-field valence groups`，隐藏按钮 tooltip 也只解释配体/d 层。`src/orbitals/mo_diagram_chemistry.cpp` 2647–2651 摘要分支把 PiFamilyOnly/Multicentre 以外的所有模式都映射为配体场。集中修复按实际模式使用普通价层或配体场措辞，补四模式及四语测试，不改变科学筛选。
- 几何独立参考：NIST CCCBDB https://cccbdb.nist.gov/expgeom2x.asp?casno=74851&charge=0 ，实验 D2h、C=C 1.339 Å、C–H 1.086 Å，与当前坐标相符；固定几何单点不能声称已通过该计算方法的优化/频率。D2h 字符表：https://www.staff.ncl.ac.uk/j.p.goss/symmetry/D2h.html （2026-09-05）。
- 证据：`OLD-009-loaded.png`、`OLD-009-structure.png`、`OLD-009-mo-normal.png`、`OLD-009-homo.png`、`OLD-009-lumo.png`、`OLD-009-diagram.png`、`OLD-009-compact-header.png`、`OLD-009-expanded.png`、`OLD-009-compact-restored.png`、`OLD-009-export-state.png`、`OLD-009-ocr.json`。OCR 字符错读已按实图和原始数据校正。当前不通过，不把标签正确等同于整体通过。

### OLD-010 二硼烷共享桥键子空间

- 8 原子、30 壳、58 AO/MO、0/singlet、8/8，producer D2H。最低 opacity=0.02 的结构图显示两个 B、六个 H，没有错误虚线；默认实际连接数据为八条 B–H 普通实线，四个 terminal 距离 1.196 Å、四个 bridging 距离 1.339 Å。部分远侧连接受当前视角遮挡；尝试拖动后图像未改变，不能声称已旋转检查所有遮挡边。`structure-rotated` 文件不作成功旋转证据。
- 实际 source3 为 Ag、-0.675153541 Ha、2e；source5 为 B1u、-0.435883078 Ha、2e，可分别选取；原始 GBasis 波函数对 D2h 操作的独立标签一致。两者 tooltip 是 `2x3c2e (shared 4c/4e source)`，confidence=64%，不是两条独立 canonical MO 各绑定一个桥氢。
- 紧凑 2 行/4e，展开 13 行（6 占据、7 虚），恢复 2 行。源码对等价桥氢联合参考空间只选择 occupied projector subspace；独立 Löwdin 诊断也显示 source3/5 的桥氢权重明显，分别约 32.4%/43.9%。虚 MO10/20 仍有较高桥氢权重，后续须按完整子空间投影确认哪些虚态应保留，不能仅凭能量或固定编号补轨道。当前只确认“占据桥键源空间”，不宣称完整成键/非键/反键活性空间已验证。
- 实验几何参考 NIST CCCBDB：https://cccbdb.nist.gov/expgeom2x.asp?casno=19287457&charge=0 ，D2h，terminal B–H 1.200 Å、bridge B–H 1.320 Å。当前固定几何 SP 与实验接近，但无同级优化/频率，不能据此证明极小结构。桥键 canonical 混合参考 Wang/Pang/Huang (2006), https://doi.org/10.1016/j.elspec.2006.01.003 ，仅作定性解释，不把异方法 MO 编号当真值。
- EXP-001 再现：实际 GUI 导出 16 行/16 MO、valence-central、Unknown/null 电子态，而屏幕紧凑 2 行；实际包 `OLD-010-original.mo.*` 已归档。展开摘要也复现 UI-001 的配体场措辞。
- 证据 `OLD-010-loaded.png`、`OLD-010-structure.png`、`OLD-010-mo-normal.png`、`OLD-010-compact-header.png`、`OLD-010-compact-complete.png`、`OLD-010-bridge-b1u.png`、`OLD-010-bridge-ag.png`、`OLD-010-expanded.png`、`OLD-010-compact-restored.png`、`OLD-010-export-state.png`、`OLD-010-ocr.json`。OCR 已按实际图/输入校正。当前不通过（共同 NUM-001/EXP-001），桥键虚空间完整性待核验。

### SYM-003 — Gaussian 线性群缩写 SG 被错误排成 S 的 g 下标

- OLD-011 一氧化碳的 producer 为 C*V，原始 `SG` 表示 Sigma，不是具有反演偶性的 S_g。界面 source7 的对称性排版却呈现 S/g 下标，原因位于 `mo_diagram_v2_part4.inc::parse_symmetry_notation`：只识别完整 sigma/pi/delta/phi 前缀，SG 落入一般字母/后缀解析。
- 独立原始基组操作检验：对分子轴旋转 60°、45°，及包含轴的镜面反射；反演不满足同元素几何置换，不添加 g/u。28 个 MO 都有一致的 Sigma+、Pi 或 Delta 标签，被测最大相对残差 3.11543e-8。source7 为 Sigma+，source5/6 和 8/9 均为 Pi。
- 待集中修复：规范 Gaussian 线性群缩写到内部表示的转换，保留 producer 原字符串和来源；盘点 SG、PI、DL 及真正的 g/u、反射正负后缀，不作全局字母替换。加入异核/同核、旋转几何、简并空间测试。

### TOPO-001 — 多个正交 π 通道被错误称作 spiro 网络

- OLD-011 CO 的 source6/8/9 tooltip 给出两个正交通道、2 原子、4e，却标记 `orthogonal spiro pi network`。二原子体系没有环，不能有螺环。
- 源码 `orbital_chemistry.cpp` 的 topology 分支把 `bundle.spiro || bundle.networks.size() > 1u` 都映射为 Spiro。独立通道数量不是螺环拓扑判据；还须审查 bundle.spiro 的上游赋值，不能仅改文字。
- 待集中修复：分开“轨道方向通道”和“原子连接拓扑”；螺结构必须有真实环且恰共享一个原子。CO、乙炔、累积双键是负例，真正螺共轭体系为正例；连通性、共轭网络与标签一起回归。
- 定义参考 IUPAC Gold Book：https://goldbook.iupac.org/terms/view/S05881 ，https://goldbook.iupac.org/terms/view/S05884 （2026-09-05）。

### OLD-011 一氧化碳

- 2 原子、12 壳、28 AO/MO、0/singlet、7/7；同源单点日志收敛、无优化/频率或稳定性证据。C–O=1.128 Å 与 NIST CCCBDB 实验几何相符：https://cccbdb.nist.gov/expgeom2x.asp?casno=630080 。不能因此宣称同级极小已验证。
- opacity=0.02 下有一条普通 C–O 连接。当前普通 ball-and-stick 渲染固定使用单圆柱，并非 Lewis 键级图；独立 Mayer 约 2.52088，不能仅凭单圆柱判为漏掉三键。
- 正常表面 opacity=0.92 / isovalue=0.03 / grid128；source7 为 -0.387651465 Ha / 2e，sigma 98%，bonding unavailable/0%，不能把 unavailable 偷换为确定非键。source6 为 -0.494592108 Ha / 2e，source8/9 为 -0.0180630141 Ha / 0e；两组 Pi 的两个成员均可分别选择，已检查的选择没有改变筛选集合。
- 紧凑和展开都为 5 行（7 个 MO），恢复仍为 5 行；相关图的上下部分分别留有截图。真实 GUI 导出为 13 行/13 MO、valence-central、Unknown/null 电子态，继续复现 EXP-001。四个本轮新导出已移到 `OLD-011-original.mo.*`，未覆盖原始输入。
- 证据 `OLD-011-loaded.png`、`OLD-011-structure.png`、`OLD-011-mo-normal.png`、`OLD-011-homo.png`、`OLD-011-pi-b.png`、`OLD-011-lumo-a.png`、`OLD-011-lumo-b.png`、compact/header/expanded/restored/export-state 图和 `OLD-011-ocr.json`。OCR 错读按实图与原数据校正。当前不通过（NUM-001、EXP-001、SYM-003、TOPO-001），普通 sigma MO 成键性质仍待独立定量判定。

### OLD-012 六羰基铬

- 13 原子、83 壳、199 AO/MO、0/singlet、54/54，producer OH；固定几何单点有同源收敛日志，没有同级优化/频率证据。最低 opacity 下看到六条 Cr–C 配位联系与六条 C–O 普通连接；Cr–C=1.916 Å，C–O=1.171 Å，独立 Mayer 分别约 1.22292/2.39852。
- 紧凑图 8 行、18 个 canonical MO（9 占据、9 虚，18e）；原始编号为 34 A1g、47–48 Eg、49–51 T1u、52–54 T2g、55–57 T1u、61–63 T2g、64 A1g、68–69 Eg。独立原生基组波函数对 Oh 旋转与反演的字符检验与这些标签一致，被测最大相对残差 4.11534e-8。字符表 https://www.staff.ncl.ac.uk/j.p.goss/symmetry/Oh.html 。
- source53 为 T2g、-0.2753064500 Ha、2e、组占据 6e；d 60.9%、配体 p 32.5%、π 99.8%、M–L overlap +0.027265。source63 为 T2g、-0.0182660266 Ha、0e；d 32.4%、配体 p 59.8%、π 100%、M–L overlap -0.055202。可分别选取，切换保持筛选集合。
- Δπ tooltip=0.2570404233 Ha、pair confidence=0.800，与上述两组平均能量差一致；不能将这一成对能量差直接称为配体场 Δo。source53 的 bonding unavailable/0% 与 pair confidence 80% 是不同证据字段，源码成键提升阈值未满足，暂不把 unavailable 本身判为错误非键；其化学分类完整性留待独立定量分析。
- 展开 25 行，恢复紧凑回到 8 行并保留 source63 选择。实际 GUI 导出仍为 valence-central、27 行/27 MO、Unknown/null 电子态，再现 EXP-001。本轮四文件已归档 `OLD-012-original.mo.*`，未改原输入。
- 另记 UI-002（Minor）：摘要 `visible levels=8/199; occupied=37` 混合最终可见行数、总 canonical MO 数和筛选前占据数。展开/紧凑仍显示相同 occupied=37，不能当成当前图占据数。集中修复需明确行/MO/占据的计数范围，不据此改电子数或筛选。
- 独立 Oh 首轮还覆盖 OLD-003/076/129/135/165/166/268/273；CrCO6 的 164–169 是能量相近的六维混合空间，不符合单一 Oh 不可约表示，不能为了填标签强行命名。其反演 trace=0，后续需拆成 T2g/T2u 子空间验证。原子病例不套用分子 Oh 分类。
- 证据：`OLD-012-loaded.png`、`structure`、`mo-normal`、`t2g-b`、`t2g-star-c`、`splitting`、compact/expanded/restored/export 图和 `OLD-012-ocr.json`；OCR 错读已按原图/原始数据核对。当前不通过（共同 NUM-001/CHEM-001/EXP-001），不把局部标签与数值吻合等同于完整科学通过。

### TOPO-002 — 双核体系的远侧配体被补成第一配位层连接

- Major；OLD-013 Re2Cl8(2-)。最低 opacity 的实图及生产连接快照均出现 17 条边：1 条 Re–Re、8 条近侧 Re–Cl、8 条额外远侧 Re–Cl。近侧距离 2.343 Å / Mayer≈0.99457，远侧 3.597 Å / Mayer≈0.17437；Re–Re=2.220 Å / Mayer≈2.97735。不能把非零远程密度耦合直接当第一配位层。
- 结构依据 Cotton/Harris (1964), Science 145, 1305，https://doi.org/10.1126/science.145.3638.1305 ，PubMed 摘要 https://pubmed.ncbi.nlm.nih.gov/17802015/ 明确描述无卤素桥支撑的双核金属键；仅查得摘要，不声称读过付费全文。
- `interaction_graph.cpp` supplementary coordination 路径使用 Mayer>=0.02、距离<=1.65 倍半径和遮挡判断，但 `covalent_adjacency` 不含已排入队列的强金属—配体连接，导致更近的 Re 不能参与遮挡排除。集中修复须建立可靠近邻层再处理补充联系，同时保留真实 μ-Cl、非对称桥、弱配位和高配位数正例；禁止仅整体提高阈值而重现 28/31 漏键。

### SYM-004 — 局部配位场标签覆盖有效的全分子 MO 标签

- Major；OLD-013 原始 D4h 标签存在，轨道浏览器声明 216 producer 标签；MO 图却把 HOMO84 B2g、LUMO85 B1u 都显示为局部 C4v 的 B2，82/83 Eu 与 87/88 Eg 都变成 E。局部五配位 C4v 近似本身不等于全分子 D4h，也不应覆盖后者。
- `mo_diagram_chemistry.cpp::recover_local_ligand_field_symmetry` 无条件写入已有 `metadata.symmetry`。须分别保存全分子标签、局部投影标签、点群/原点/坐标系/来源/残差，筛选及配对明确使用哪个作用域。
- 独立 `audit_symmetry.py --group D4H` 从原生 Gaussian 基组与系数计算 C4z/C2z/C2x/C2xy/反演；上述前线表示全部通过，残差约 1e-14 至 5.45e-11。字符表 https://www.staff.ncl.ac.uk/j.p.goss/symmetry/D4h.html 。同类计算已覆盖 OLD-077/098/127/248，尚不等于这些案例 UI 通过。

### CHEM-002 — 金属—金属通道被全分子汇总掩盖，且错误套用 donor 配对

- Major；OLD-013 底层 source84/85 的 Re–Re pair channel 已为纯 δ，但 MO 图仅突出全分子 π 汇总。独立 GBasis overlap 和对称 Löwdin 投影得到：84/85 的两个 Re 总权重 0.685896/0.778517，均为相对 Re–Re 轴的 d|m|=2；每单位占据 pair overlap 为 +0.0468588/-0.0584836。这支持 M–M δ/δ*，同时允许 M–Cl π 尾部共存，不能只保留一个无作用域 family。
- source82/83 与 87/88 的 Re dπ 权重约 0.406974/0.825742，Re–Re overlap 为 +0.0891939/-0.660976。UI 却把全分子 Eu 与 Eg 配成 `pi-donor splitting`，Δπ=0.1950575785 Ha，confidence=0.696。能量差算对不代表 donor 机制正确；不同全分子 irreps 的 M–M 成键/反键对不能只因局部同属 E 就视作同对称的 M–L 混合对。
- 集中修复：区分 M–M、M–L、配体内部相互作用；同一 MO 可以有多个带原子/轴/证据的通道和成键描述。M–L donor/acceptor 配对须满足其全分子/适用子空间对称性及局域成分证据，M–M 对另行标注，禁止用固定分子名或 MO 号补 δ。
- 独立结果 `OLD-013-pair-reference.json`；生产 pair 诊断 `OLD-013-cov-metal-pair.json`。Löwdin population 和 Mulliken overlap 都依赖基组/划分，不冒充唯一可观测键级；本例 S 条件数约 38862.65，虚轨道过大 Mulliken 数值另作方法限制，不能凭数值大就判坏波函数。

### OLD-013 双核铼首轮 UI 与导出

- 10 原子、108 壳、216 AO/MO、-2/singlet、84/84（含 ECP 的有效电子计数）；有同源收敛单点日志，无同级优化/频率或稳定性证明。
- source84/85、82/83、87/88 均可单独选择；紧凑 18 行/21 MO/28e，展开 40 行，恢复 18 行并保留 source88。实际导出 25 行/25 MO、valence-central、Unknown/null 电子态，再现 EXP-001；本轮四文件已移入 `OLD-013-original.mo.*`，没有覆盖输入。
- 证据 `OLD-013-loaded.png`、`structure`、`mo-normal`、`homo`、`lumo`、`eu-a/b`、`eg-a/b`、compact/header/expanded/restored/export 图及 `OLD-013-ocr.json`。OCR 对 84/88、能量和 g/u 存在错读，结论已按实际截图和原始数据校正。当前不通过；更深的 Re 半芯层筛选仍待独立径向/价层判据，不能单凭深能量删轨道。

### UI-003 — 开壳层合并行的能量、成员选择与组成信息来源混用

- Major；OLD-014 的三个 π 空间行合并 α11/β11、α12/β12、α13/β13，但分别仍使用 α 能量 -0.337909924、-0.214322853、+0.0431125195 Ha。对应 β 能量为 -0.302741300、-0.066138023、+0.0726226343 Ha，不是同能简并。
- 从原始轨道列表搜索 79 并选择 β12 后，再点击紧凑图中间行，tooltip 正确显示 source79/internal78、Beta、0e、-0.066138023 Ha；水平线却仍在 -0.214322853 Ha。tooltip 的 ligand-p=92.5% 也沿用 α12 的行 metadata，而独立的生产 β12 channel/manifold 为约 88.9%。不得让单个 β MO 的能量/占据与未注明来源的 α 组成拼成一份说明。
- `mo_diagram_chemistry.cpp::merge_spin_counterparts` 保留 α 的 layout_energy/metadata；`orbital_ui_v2.cpp` 的 member_points 只有当 counterpart 已在别处选中时才返回 β，否则返回 α，图内没有分别进入两套自旋的明确操作。须设计紧凑的双自旋选择与能量呈现，不可用简单全展开破坏整洁，也不能把 α/β 当作同一能量。
- HOMO/LUMO 快捷按钮按当前选中的自旋导航，源码明确如此；本次 α 选中时 LUMO 跳到 α13，不另记为排序错误。修复时可明确该自旋上下文，避免误解为跨自旋最低空轨道。
- 证据 `OLD-014-beta12-browser.png`、`OLD-014-beta12-compact-tooltip.png`，与先前 `OLD-014-somo-alpha.png` 配对；原始索引、浮点值以 FCHK 和实图为准，OCR 存在 79/712 等错读。

### OLD-014 烯丙基自由基首轮 UI 与独立分析

- C3H5，8 原子、33 壳、67 AO、134 自旋 MO，0/doublet、12/11；同源 NoSymm UDFT 单点收敛，S² before/after=0.7929/0.7503，没有同级优化、频率或稳定性证据。最低 opacity 下 7 条普通实线连接正确可见；两条 C–C=1.390 Å、Mayer≈1.51293，五条 C–H≈1.080–1.090 Å、Mayer≈0.959–0.966。未复现含氢虚线错误。
- 对称性缺失扩展 SYM-002：COV 检测 C2v，但 134 个 α/β 标签全不可用。独立原生基组/系数的 C2、两镜面操作可恢复全部 134 个标签；前线 α/β 11、12、13 均依次为 B2、A2、B2，前线残差不超过 5.82e-14。本参考以 C2 为 z、分子面为 xz；B1/B2 必须随坐标约定解释，不以标签交换误报。字符表 https://www.staff.ncl.ac.uk/j.p.goss/symmetry/C2v.html 。平面 C2v 参考已额外处理同类候选，但不适用非平面/无法确定轴的结果不能判成输入坏文件。
- 紧凑 3 行包含 α/β 各三个 π MO、合计 3e、一条垂直分子面的通道。展开 16 行，恢复 3 行，α 成员均可选择。β 单独选择及图中能量问题见 UI-003。正常表面继续受 NUM-001 影响，不能凭现有颜色形状断言波函数错误。
- SOMO 的 anti/95% 保留为需修正作用域的化学解释问题：独立相邻 C1–C2 overlap 为 α12 +0.007524963、β12 +0.008547234，而端点 C1–C3 为 -0.08386326/-0.09966549。生产实现因中心原子权重门限把相邻 pair 记为 NotApplicable，端点反相 pair 主导全局 anti。与相邻键近似非键的说法不矛盾；应分别说明近邻共轭键与远端相位，不能直接把负 overlap 数值删掉，也不能只凭最小 Hückel 图强改成绝对非键。证据 `OLD-014-pair-reference.json`、`OLD-014-terminal-pair-reference.json`；与 CHEM-002 一同处理带作用域的角色描述。
- 实际 GUI 导出为 25 行/25 MO、valence-central、Unknown/null 电子态，继续复现 EXP-001；本轮四文件已移动到 `OLD-014-original.mo.*`，未覆盖原输入。structure、mo-normal、三个 α tooltip、β12、compact/header/expanded/restored/export 截图及 `OLD-014-ocr.json` 已保留。当前不通过，没有完整科学通过结论。

## 严格回归报警（尚不能直接判定软件或断言哪方错误）

| ID | 原库案例 | 报警 | 下一步独立判定 |
|---|---|---|---|
| REG-001 | OLD-003 TiF6(2-) | compact T1/T1g exclusion damaged the metal framework | 比较当前筛选的 σ 框架、d 子空间与原始权重，检查断言是否过时 |
| REG-002 | OLD-020 cp_anion_pi_vmd | 单成员行 declared degeneracy=2 | 检查筛选拆散的简并组、producer/derived 对称性及公共 builder |
| REG-003 | OLD-238 MoCN8(4-) | square-antiprismatic CN8 open-shell regression | 独立 α/β 子空间、局域 D4d 与全分子点群区分 |
| REG-004 | OLD-240 HfH2O10(4+) | pentagonal-antiprismatic CN10 regression | 第一配位层、D5d 局域几何、紧凑选择覆盖 |

## 输入质量与证据边界

- OLD-002 MO11（raw index 10，-0.398842813 Ha，2e，Sigma_g）UI 标为 sigma/antibonding、89%，三中心 3c4e 93%。暂记化学归属待核验，不因教材最小三轨道模型常称中间成员“非键”就直接判软件错；需看真实 Xe 5s 混合、M–L overlap 与分类公式。MO6（raw index 5，-0.548762397 Ha）为 Sigma_u、bonding 97%。原文件为 SDD Cartesian，抽测 CUDA 网格通过。

- 273/273 文件独立解析、MO 重建 density 与 producer 检查已完成。
- Multiwfn 2026.8.21 对全部 273 的每个 α/β 轨道归一化及 Tr(PS) 通过；最大 norm 偏差 1.393e-7。UHF 必须使用能输出两个自旋块的 GTF-form 测试，不能把只打印 alpha 的测试算成 beta 通过。
- 另外使用 IOData 1.0.1 / GBasis 0.1.0 原生解析 AO 积分已完成全部 273 的完整同自旋 CᵀSC 检查，273 通过；alpha/beta 分开检验，不要求不同自旋空间互相正交。逐例误差、S 条件数和 Tr(PS) 在 `overlap_cases.jsonl`。COV 的反推 S 还需与这个独立 S 比较；这不取代已发现失败的 CUDA 网格检查。
- OLD-019 `cp_anion_pi.fch` 缺明确同源日志。其余 272 已核对；仅 23 有优化完成和频率证据。单点/扫描不自动判错，不宣称已证明极小或正确基态。
- OCR 原始识别包含撇号/百分比误读；本台账轨道标签数值以截图人工回看及输入核对为准。最低 opacity 结构证据与正常 opacity MO 证据分开保存。
- COV 与独立 S/Mayer 的全库对照已完成：272 数值符合当前容差，OLD-020 的截断 MO 例外见 NUM-002；最大绝对 S 差 6.15913e-6、Mayer 差 1.33588e-6。分布不变量发现 260 文件、23,850 个 MO 成键分布异常（CHEM-001）；正确允许 NotApplicable 的零分布后，channel/manifold 无同类报警。原始结果以 schema_version=2 为准，不能把 reference_status 当 COV 整体通过。
- Gaussian 新计算 0；OLD-001 至 OLD-014 已收集首轮主要 UI/导出证据，尚没有完整科学通过案例；不得把本台账称为全量完成报告。

## 2026-09-05 原生验证首例 OLD-015 萘（新增执行记录）

- 分支 `test/fchk-validation-native`，本地提交 `cc08ecdc551809d4f5d61b4a55a41e9afca24da8`；生产科学算法未修。输入为 `15_naphthalene_D2h_def2SVP.fch`，18 原子、84 壳、180 AO/MO，0/singlet，34/34。同源 PBE1PBE/def2SVP 单点日志正常结束；无优化、频率、稳定性证明，不把输入自动判坏或宣称稳定基态。
- 首轮执行及图像复核完成，整体不通过。报告和完整证据位于 `F:\Codex\2026-09-05\branch-15\outputs\OLD-015-验证报告.md` 及 `OLD-015-naphthalene\`。自动耗时 35.189 秒（界面/GPU 29.498 秒）；包括复核、分析、报告准备的端到端耗时 946.660 秒；首次设施开发及编译另计。
- 独立 CᵀSC 最大归一化误差 1.48e-8、off-diagonal 3.66e-8，Tr(PS)=67.9999999766；COV S 最大差 1.78e-7，33 对 Mayer 最大差 1.98e-8。结构 19 条普通共价连接（11 C–C、8 C–H），两个稠合六元环；蓝色 π/芳香辅助虚线另行显示。
- 紧凑 π MO 为 27/30/32/33/34/35/36/37/41/47（从 1 编号），5 占据+5 虚，共 10e；独立 D2h 与 10 个标签一致。全 180 个中独立参考恢复 178 个；近核心 MO5/6 的近能子空间留为参考边界。10 个成员均有实际选择和图像，10→25→10 恢复正确。673 动作、0 采集失败，4574 帧；已区分 191 个正常延迟更新帧，当帧场景/UI 错配 0。
- **NUM-001 继续复现于真实渲染器。** 每个 MO 8179 个确定性去重纹理样本，180/180 不通过。HOMO34 NRMS=0.0465895，LUMO35=0.0738169，最差 MO180=2.4163053；整体正负相位对齐仍失败。独立参考副本的纯壳奇数绝对磁量子数相位实验使全部 180 MO 的 NRMS 降至 5.48e-6 以下；实验没有修改生产值，不算修复通过。
- **CHEM-001 继续复现。** 172 个 MO 成键分布非法；HOMO34 四项为 [0, 0.9253585401, 0.0746414599, 15.2898047512]，UND=1528.98%、总和=1628.98%。界面分类置信度 96% 是另一字段，不与 92.54% 反键分布混为同值。
- **EXP-001 继续复现。** GUI 为 10 行紧凑 π 空间，实际 PNG/SVG 和 JSON 布局为 25 行 valence-central，成员不一致；电子态 Unknown/null。JSON/CSV 的 180 轨道记录及 SVG 布局元数据自身一致不能抵消导出与 GUI 不一致。
- **CHEM-002 作用域补充。** HOMO34 UI 为全局 antibonding / 96%。独立相邻 C–C Mulliken overlap 每单位占据合计 +0.217597669，8 正/3 负（范围 -0.0446602 至 +0.0747954）；该数值依赖划分，不能冒充唯一可观测键级或把所有占据 π 都强行叫成键。继续要求局部键/远程相位/全局汇总分别给出作用域。原始证据 `pi-pair-reference.json`、`native/pi-034-tooltip.png`。

### UI-004 — 长轨道提示框越出窗口底部，化学说明无法完整阅读

- Minor；OLD-015，在 2100×1250 帧缓冲、UI 缩放 1.5 下复现。MO34 提示框从接近顶部铺至底边之外，参与原子与后续化学解释被裁切；其他 π tooltip 也有长内容风险。证据 `F:\Codex\2026-09-05\branch-15\outputs\OLD-015-naphthalene\native\pi-034-tooltip.png`，同目录保留各成员图与 UI 数据。
- 当前 `draw_level_tooltip` 把完整内容堆入无界高度 tooltip；本轮只登记，未修生产界面。后续应约束可用高度并提供可访问的滚动或独立详情，同时保留轨道身份与关键数据；不能以后台 JSON 存在替代用户可读性。
- 图像复核不因此声称裁切部分已经可读；四语言仅验证切换与已见内容，完整术语未认证。本例已完成首轮但不通过，后续 OLD-016 尚未开始。

## 2026-09-05 全库纠错分析完成检查点（历史版本，裁决见后附源码设计复审）

- 用户本轮已授权分析取证内容、逐条登记、分类及提出通解规划。本检查点取代此前“等待纠错分析”的状态；当前仅完成分析与规划，没有修生产算法或启动新 Gaussian。
- 冻结验证分支 `test/fchk-validation-native`，提交 `2284ab4b048c526b468f026ddf8fadaf296df10f`，工作树干净。273文件、31,871实际MO纹理、951,491帧、11,378PNG均纳入对应分析；真实结构联系表273例已概览，自动检查与人工覆盖边界分开记录。
- 新的全MO采样结论：NUM-001涉及252文件/26,229MO；CHEM-001涉及260文件/23,850MO；EXP-001全273例电子态丢失，272例成员集不同、90例模式不同。此前247文件网格报警是旧代表轨道抽样口径，本轮全MO口径取代影响范围，不能混写为同一次检查。
- 对称性独立可分类28,902MO；641个错误derived标签（520偶数Dnh、121群上下文混用），13,819个非平凡可计算空标签；2,969MO参考未定保留。合法别名/轴约定/C1平凡缺省不计真错误。
- 21个台账项、67,360条发生记录：17软件问题、1功能覆盖缺口、1取证视野问题、2测试规则问题。每条记录关联输入身份、原MO/spin或原子对/图像及具体证据；同一案例可有多项错误。
- 四项REG报警已判：REG-001为T1G/T1g字面比较误报（QA-001）；REG-002为18行单成员仍声明二重简并的真问题（UI-007）；REG-003/004分别为12>11、10>8的固定行数政策（QA-002），并不能证明其化学正确或错误。其他独立失败仍保留。
- 新增细分：SYM-005群上下文不一致；UI-005固定视口遮挡；UI-006近似局域组误称严格简并；UI-007拆组元数据残留。4个初次compact集合报警经counterpart映射解除；OLD-107真螺环列为保护负例。
- 最新完整报告/可搜索逐例入口：`F:\Codex\2026-09-05\branch-15\outputs\fchk-analysis-20260905\index.html`；主报告`F:\Codex\2026-09-05\branch-15\outputs\fchk-analysis-20260905\全库取证分析报告.md`；根因台账`F:\Codex\2026-09-05\branch-15\outputs\fchk-analysis-20260905\问题分类与根因台账.md`；通解规划`F:\Codex\2026-09-05\branch-15\outputs\fchk-analysis-20260905\错误分类与通解规划.md`；逐条记录`F:\Codex\2026-09-05\branch-15\outputs\fchk-analysis-20260905\逐条问题记录.jsonl`。
- 执行顺序：AO约定/S → 分布/图拓扑/化学作用域 → 群表示及global/local来源 → 自旋/近似组/筛选/导出统一数据 → 全273例独立参考与实际UI/普通OFF构建联合复验。禁止分子名、固定MO编号、原子序号特判或删轨道凑行数。
- 尚未修复：全部273例科学判定仍不通过（每例至少EXP-001）。独立活性空间、二硼烷虚空间、未定irrep、完整体素和遮挡/语言细节按H/G类保留，不能用取证完整宣称全库科学通过。历史输入、旧证据、原用户工作树修改及后续至少50新分子目标均保留。

| 问题项 | 影响文件 | 逐条记录 |
|---|---:|---:|
| NUM-001 — 纯球谐 AO 相位约定不一致 | 252 | 26,229 |
| NUM-002 — 矩形 MO 系数块缺失 overlap 和 Mayer | 1 | 1 |
| CHEM-001 — 成键分布累计器导致百分比超过 100% | 260 | 23,850 |
| EXP-001 — 实际导出与当前界面构建来源分裂 | 273 | 273 |
| SYM-001 — 偶数阶 Dnh 派生标签使用撇号并缺 B 表示 | 7 | 520 |
| SYM-002 — 独立可恢复的非平凡对称性标签缺失 | 160 | 13,819 |
| SYM-003 — Gaussian 线性群缩写没有正确规范化排版 | 2 | 6 |
| SYM-004 — 局域配体场标签覆盖 producer 全分子标签 | 1 | 19 |
| SYM-005 — 派生 irrep 与展示的全分子点群不一致 | 3 | 121 |
| TOPO-001 — 多个正交 π 通道误称螺环 | 26 | 81 |
| TOPO-002 — 双金属远侧配体被补成第一配位层 | 1 | 8 |
| CHEM-002 — 全分子、原子对与局域子空间的成键语义混用 | 4 | 9 |
| UI-001 — 非金属价层摘要被称为配体场 | 133 | 133 |
| UI-002 — 摘要占据/虚轨道计数与当前可见集合不同 | 107 | 107 |
| UI-003 — α/β 合并行混用能量和组成 | 31 | 624 |
| UI-004 — 长提示框底部裁切 | 143 | 1,448 |
| UI-005 — 结构取证视野被左侧面板遮挡 | 28 | 29 |
| UI-006 — 近似局域分组合并被当成精确简并 | 12 | 62 |
| UI-007 — 拆分成单成员后仍保留二重简并元数据 | 1 | 18 |
| QA-001 — TiF6 回归对标签大小写产生误报 | 1 | 1 |
| QA-002 — 固定行数上限被误作化学错误判据 | 2 | 2 |


## 2026-09-05 源码设计复审与台账裁决检查点（当前最新）

- 用户要求：读上传的源码规则分析与引用任务“本地继续 COV”，依据真实设计意图区分合理近似、需要调整的政策、明显错误，修订既有台账。本阶段完成分析与规划，没有实施生产修复或新Gaussian计算。
- 冻结验证分支：`test/fchk-validation-native`；HEAD `2284ab4b048c526b468f026ddf8fadaf296df10f`。源码保持干净。附件引用的15个文件与当前基线归一化换行后相同，不能把CRLF/LF造成的不同哈希说成源码变更。
- 共复审55组规则：27组保留设计，16组调整设计，7组涉及确定错误，5组属于能力/证据边界。规则组不是缺陷数。
- 旧21项重新裁决＋新增3项＝24项：10项确定缺陷、9项设计/表达修订、2项能力缺口、1项证据限制、2项测试问题。各项可能关联多条相关记录，不能相加为独立化学错误。
- 67,360条旧记录全部逐条迁移裁决，原文与原证据未删。Zn近似非键1条排除错误；萘/烯丙基2条撤销确定全局成键判错。624条自旋代表线差异和62个近似局域组改为呈现/作用域问题；不是624条原始β能量错误或62组源能量被改坏。
- NUM-002改能力缺口；SYM-002保持能力缺口。TOPO-002撤销“远侧8条进入第一壳/CN被多算”，实际主壳4Cl＋另一Re、CN5，保留次级联系分层修订。UI-007明确源二重关系和单行成员，不判源简并错误。QA-001/002排除产品化学报警但保留测试修订。
- 新增POL-001：OLD-004 E/T2晶场间隔被放入π分裂字段，1条实际配对观察；保留用户许可的近似非键。新增POL-002：全壳SigmaOnly硬否决遇P配体先验过宽，有源码/原始研究反例，本库P配体场样本0，漏检量未确定。新增UI-008：OLD-020分析不可用却显示组成/通道默认零，23条实际tooltip观察。这24条来自已冻结取证，不是新运行。
- SYM-004的19个metadata替换再次核对；真实回写点为`ligand_field.cpp:908`，不是最终仅补空值的循环；源FCHK/wavefunction全局标签未被覆盖。CHEM-002另核实局域M–L聚合可包含第一壳的另一个金属，需要明示中心邻域/MM/ML作用域。
- NUM-001（252文件/26,229MO）、CHEM-001（260文件/23,850MO）、EXP-001（273文件）等确定缺陷仍未修复。无新增完整通过案例。

| ID | 本轮操作 | 当前类型 | 严重程度 | 当前标题 |
|---|---|---|---|---|
| NUM-001 | 保留 | 确定缺陷 | Critical | 纯球谐 AO 约定错配 |
| NUM-002 | 降级 | 能力缺口 | Capability | 矩形 C 的 S 恢复能力不足 |
| CHEM-001 | 保留 | 确定缺陷 | Moderate | 成键分布累计器初始化错误 |
| EXP-001 | 保留 | 确定缺陷 | Major | 实际导出未使用当前界面数据 |
| SYM-001 | 保留 | 确定缺陷 | Major | Dnh 派生标签语法/分类错误 |
| SYM-002 | 边界保留 | 能力缺口 | Capability | 非平凡全局 irrep 派生能力缺口 |
| SYM-003 | 保留 | 确定缺陷 | Moderate | Gaussian 线性标签排版错误 |
| SYM-004 | 保留并纠正根因 | 确定缺陷 | Major | 局域标签占用唯一 Symmetry 字段 |
| SYM-005 | 保留 | 确定缺陷 | Major | 派生标签与全局点群字段未绑定 |
| TOPO-001 | 保留 | 确定缺陷 | Major | 多方向 π 网络误称螺环 |
| TOPO-002 | 改写并降级 | 设计/表达需调整 | Moderate | 远侧配位联系缺少次级联系显示层 |
| CHEM-002 | 部分排除并改写 | 设计/表达需调整 | Moderate | 多中心/多键轴解释缺少明确作用域 |
| UI-001 | 保留 | 设计/表达需调整 | Minor | 普通价层模式使用配体场措辞 |
| UI-002 | 改写并降级 | 设计/表达需调整 | Minor | 摘要混合行数和不同集合的 MO 数 |
| UI-003 | 改写并降级 | 设计/表达需调整 | Moderate | 自旋代表行与单 MO 详情来源未分清 |
| UI-004 | 保留 | 确定缺陷 | Minor | 长 tooltip 不能完整阅读 |
| UI-005 | 边界保留 | 证据限制 | Evidence | 固定取证视角的遮挡限制 |
| UI-006 | 改写并降级 | 设计/表达需调整 | Moderate | 近似局域组缺少近似标注及能量范围 |
| UI-007 | 改写并降级 | 设计/表达需调整 | Minor | 单成员行与源二重组关系表述不清 |
| QA-001 | 排除产品化学报警 | 仅测试问题 | Test | TiF6 标签大小写测试误报 |
| QA-002 | 排除产品化学报警 | 仅测试问题 | Test | 固定行数测试误报 |
| POL-001 | 新增 | 确定缺陷 | Moderate | 弱 d 壳晶场间隔被显示为 π 分裂 |
| POL-002 | 新增 | 设计/表达需调整 | Major | 饱和 P 先验在全壳 SigmaOnly 时硬否决 π 配对 |
| UI-008 | 新增 | 设计/表达需调整 | Minor | 分析不可用时组成和通道仍显示默认零 |

后续通解顺序：底层AO与分布不变量 → 全局/局域、MO/组、有效0/N/A的数据约定 → π配对与晶场类型/先验权限 → 拓扑与次级联系显示 → 实际界面/导出统一 → 冻结273例及异类独立输入复验。保留紧凑、π-only、局域近似、真实简并、自旋代表行和可访问隐藏成员；不能为满足固定行数或无N/A而硬填。

当前主报告：`F:\Codex\2026-09-05\branch-15\outputs\fchk-rule-review-20260905\源码设计复审报告.md`。
逐项规则：`F:\Codex\2026-09-05\branch-15\outputs\fchk-rule-review-20260905\源码规则逐项裁决.md`。
修订台账：`F:\Codex\2026-09-05\branch-15\outputs\fchk-rule-review-20260905\修订后的问题台账.md`。
通解规划：`F:\Codex\2026-09-05\branch-15\outputs\fchk-rule-review-20260905\修订后的通解规划.md`。
逐条裁决：`F:\Codex\2026-09-05\branch-15\outputs\fchk-rule-review-20260905\逐条裁决.jsonl`。
可搜索入口：`F:\Codex\2026-09-05\branch-15\outputs\fchk-rule-review-20260905\index.html`。


### SYM-006 — Gaussian 空点群字段误读为列名 NOp

- 裁决：新增，确定缺陷，Major，尚未修复。发现时间 2026-09-06T08:26:43+08:00；已有源码问题编号保持不变。
- 全 273 例独立核骨架采集 NS-002 结束后交叉排查发现，受影响案例为 OLD-070。原日志出现 `Full point group                         NOp   0`，点群字段为空；既有真实生产取证却报告 `point_group = NOp`。
- 根因：`src/parser/gaussian_log.cpp:160` 的 `token_after` 读取空字段之后的下一列名；随后 243–253 行把它标记为 Producer。`src/analysis/symmetry.cpp:524` 遇到 Producer 即返回，进一步阻止有效的几何推导。已核对这两份生产源码相对于原取证提交 `2284ab4` 没有变化。
- 通解：按字段边界和合法语法解析原始点群，空字段、缺项、列名均明确标为不可用；保留原始文本、producer detected/used 及独立几何来源。不能为了过滤 NOp 把允许点群硬限制到当前 20 项目录。缺失/无效的 producer 字段允许几何推导，但必须保留近线性、畸变和闭合残差。
- 边界：该例核坐标在不同明确容差下呈现近线性差异，不能由这个解析错误推断输入物理错误，也不强制指定一个精确点群。回归应覆盖空字段、正常有限群、D*H/C*V 别名、跨 Link1 元数据及完全缺项，并对全正式集合复验。
- 机读证据：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\crosscheck-additions\NS-002\SYM-006.json`；完整核骨架交叉核对：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\nuclear-symmetry-reference\NS-002\post-batch-review-v2\review.json`。无既有 COV 问题在此关闭，正式完整通过数仍为 0。


### AN-001 全 MO/自旋算符作用范围补充

- 原库 273 例全部取证和诊断终止后完成本次逐例裁决，分析身份 `d50b6652ce939709a7272ca15db3b95ff4f93b91353084013617c3adc7da09ca`。12 例在核对称操作/有限探针下出现一粒子密度算符范围差异；不把它们记成 12 个 COV 缺陷。
- 排除“核骨架点群必须等于每个单自旋波函数/单 MO 标签”的机械检查。OLD-028 的 C₂ 操作交换 α/β 算符，不能仅凭核几何近 C₂v 就将日志 Cs 判错。12 个原子案例的 producer 有限 OH 子群与完整核骨架 Kh 也不能仅按字符串不同判错。
- SYM-002/004/005 的范围、推导和混合状态要求继续保留；SYM-006 的空字段 NOp 解析错误也仍是确定、未修复问题。没有关闭已确认的产品缺陷。
- 新增参考审阅覆盖问题 REFQ-003：自动候选标记 false 不能排除电子态竞争。已取得 OLD-052/098/101 的限制性反例；问题属于参考规划，不是生产 COV 缺陷。等 REF-001 全批次终止后统一审阅并冻结适用的补充作业；当前冻结输入、电子态、算法及判据未改变。
- 新证据与通解补充：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\crosscheck-additions\AN-001\作用范围裁决与通解补充.md`；完整逐例机读裁决 `F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\crosscheck-additions\AN-001\adjudication-v1.json`。密度度量为 Hilbert–Schmidt 一粒子算符范数，未宣称实空间密度 L²、完整连续群或最终物理状态通过。

- AN-001 当前机读裁决更新为 `F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\crosscheck-additions\AN-001\adjudication-v2.json`，文字报告为 `F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\crosscheck-additions\AN-001\作用范围裁决与通解补充-v2.md`。EVID-SYM-002 修正报告器的原子类型枚举匹配，使 12 个原子条目与原文“OH/Kh 范围差异”的裁决一致；原始数值、生产结果及既有产品问题状态未改变，v1 保留。


### EVID-RUN-001：验证运行器启动阶段计时缺口（2026-09-06T11:50:41+08:00，保留，诱因待定位）

SD-001 全 273 例取证结束后发现：OLD-091/092 进程树墙钟各约 500 秒，CPU 各不足 0.6 秒；恢复子进程前的启动记录比监督起始时间晚约 499 秒。两个作业最终正常完成，不能记为分子错误，也不能将等待计作密度计算耗时。通用方案为细分单调时钟阶段、对照验证辅助程序控制台模式和并发句柄继承、区分启动与执行期限。当前证据不足以认定控制台或句柄继承就是诱因。详见 [运行器计时证据](F:/Codex/2026-09-05/branch-15/outputs/cov-complete-validation-20260906/runtime-incidents/launch-latency-SD-001/启动阶段计时缺口.md)。冻结 REF-001 未修改，既有 COV 科学问题未关闭。


### EVID-RUN-002：启动阶段时钟分辨率不足（2026-09-06T12:54:43+08:00，保留，修正待复验）

LP-001 的阶段计时使用本机分辨率 15.625 毫秒的 monotonic/GetTickCount64，短操作出现零 tick，不能作为零耗时或微秒级结果。12 次协议控制结果保留，短阶段性能比较暂不接受；通解为性能计数器整数时间、显式时钟元数据及回调边界核对，随后整批复验。此问题不解释原约 499 秒等待，也不影响 FCHK 数值裁决。见 [分辨率更正](F:/Codex/2026-09-05/branch-15/outputs/cov-complete-validation-20260906/runtime-incidents/launch-phase-controls/LP-001/post-batch-review-v2/启动阶段整批对照-分辨率更正.md)。


### LP-002 高精度阶段时钟修正验收（2026-09-06T13:32:34+08:00）

EVID-RUN-002 在源码 `d802885a15eda92b467dfda9fa03a6137a26f3b8` 的阶段遥测范围内已修复并验证：3 项真实进程测试和 12 次固定对照全部通过，保留整数计数器差值、回调边界和实际时钟元数据。旧 LP-001 粗粒度记录不改写。EVID-RUN-001 仍保留，本次协调器 Popen 另测得 898.662420 秒等待，诱因未确定。未改 REF-001 冻结运行器和科学输入，不关闭 COV 科学问题。见 [验收证据](F:/Codex/2026-09-05/branch-15/outputs/cov-complete-validation-20260906/runtime-incidents/launch-phase-controls/LP-002/高精度计时修正验收.md)。


### NUM-003 — producer overlap 免检并报告未计算的零残差（2026-09-06T14:38:32+08:00，新增、未修复）

未修改的生产 enrich_fchk_overlap_from_file 在全部四项组件反例中都接受矩阵并将 Gram 最大残差记录为 0。有限矛盾矩阵实测残差 0.280637404，不定矩阵实测 1.320159351、最小特征值 −0.2，NaN 矩阵也被接受。来源优先规则保留，免检和假零值需要通解修复：从基组独立计算 S，同时验证原矩阵、全部自旋 Gram 和状态。原库 273 例的实际受影响范围尚未据此确认。见 [整批证据](F:/Codex/2026-09-05/branch-15/outputs/cov-complete-validation-20260906/regression-fixtures/producer-overlap-v1/post-batch-review-v1/可选重叠矩阵反例裁决.md)。

### EVID-RUN-003 — 编译辅助进程延长资源 Job 生命周期（2026-09-06T14:38:32+08:00，新增、通解待实现）

本次编译主进程已正常退出，独立创建的 vctip 仍存活。核对精确身份后仅清理该辅助进程，同一冻结编译作业随即结束并完成四项反例采集。记录保留此干预，不将全部编译墙钟都归因于它。后续为有限编译驱动定义受控收尾契约，Gaussian 仍以整棵进程树结束为准。EVID-RUN-001 另在本次执行中实测 CreateProcessW 阶段等待 45.7067343 秒，内部诱因继续未定。


## 2026-09-06 EVID-RUN-003 通用编译收尾验收

- 源码：验证分支提交 `fb3a24a6b6337093f20eb0f7d4bdec944a9ad372`。此提交修改通用进程监督设施及真实进程控制测试，未修改生产 COV 科学算法；正在运行的 REF-001 冻结副本保持原样。
- EVID-RUN-003 更新为 **在显式有限编译任务范围内已修复**。调用方声明可在主程序退出后结束的辅助程序完整路径；监督器通过主程序句柄确认退出，对实际 Job 成员持有句柄，核对完整映像路径、任务归属并记录创建时间。未知成员不被该收尾机制结束，整个 Job 未清空前不释放配额。Gaussian 等科学启动器继续默认等待整个任务。
- 第一轮 8 项中 3 项失败，原因是“所有残留进程都必须是声明辅助程序”的限制被附着的控制台进程阻挡。完整失败结果保留于 `runtime-incidents/finite-driver-preflight-v1/post-batch-review.json`。对应预备编译轮 producer-overlap-v2 没有启动，未计作已运行。
- 第二轮全部 9 项真实进程检查通过，0 跳过，耗时 105.722 秒。验证默认等待实际工作子进程、成功/失败和 259 退出码保留、同名不同路径保护、任务外同名进程保护、未声明工作进程继续完成，以及原有高精度计时/句柄/超时行为。没有放宽期限。
- 真实 MSVC 编译退出 0、未超时，自动收尾 1 个经核验的 vctip 辅助程序。编译阶段总耗时 89.984 秒，累计 CPU 6.922 秒；从观察到主程序退出至 Job 清空约 0.501 秒。实际控制采集耗时 0.534 秒。EVID-RUN-001 的系统启动延迟根因仍未确定，不能把全部耗时都归于辅助程序。
- 同一批四个生产 overlap 控制全部采集后再比较，所有观察字段与旧轮一致。**NUM-003 仍为已确认、未修复**；没有新增正式分子通过，也没有新增外部分子。
- 可读验收：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\regression-fixtures\producer-overlap-v3\post-batch-review-v1\编译辅助程序通用收尾验收.md`；机读证据：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\regression-fixtures\producer-overlap-v3\post-batch-review-v1\review.json`。
- 新增实施前契约 `general-fix-contracts/decomposition-v1`：细化 CHEM-001 的零累加器、状态/范围传播、参考投影守恒以及权重缩放的适用前提；明确当前绝对原子权重并非 Löwdin 布居。契约尚未实施生产修复，不关闭 CHEM-001、CHEM-002、UI-008。


## 2026-09-06 GI-001 原库结构取证及保守去重底账

- 验证分支提交 `a94135734d29edb0e4b3d19e5cc1e7572755d7c2`，增加独立的核连接图、原库批量采集和保守新增计数检查；三份受测源码与 GI-001 冻结副本逐字节一致。9 项控制检查通过，无跳过。生产 COV 科学算法、REF-001 冻结输入和运行副本均未修改。
- 原始 273 案例全部进入终止状态：273 份结构证据已采集、0 采集错误，采集主体耗时 4.190881 秒。核对原始 FCHK 哈希、核数、元素及坐标，保留全部显式氢。每例一组固定旋转/平移及确定性原子重排，在三个尺度共 819 项连接图抽查通过；这不构成 MO 或渲染不变量验收。
- 全批结束后比较 379 对同组成案例：255 对逐尺度同构、56 对逐尺度不同、68 对关系随尺度变化。25 个案例对半径尺度敏感，0 个处于 1e-6 Å 数值边界带。三个尺度得到 164/166/167 个距离图等价类，**不能当成不同化学分子数**。
- 14 组预设配对中 13 组在全部尺度符合输入意图。OLD-098/099 在 1.35 尺度的差异已查明：正方形 C—C 对角线 2.0364675302 Å 被 2.052 Å 截断包含，矩形对角线 2.0565018841 Å 则不包含；原始 13/14 结果保留，不调阈值。这两个几何不凭此差异获得新增分子计数，亦不据此新建 COV 错误。
- 6 组环连接块控制均符合预设意图。OLD-107 的两个循环块仅共享一个原子，支持保留其作为真实螺环的结构正例；本项不关闭 TOPO-001，也不替代轨道通道与实际显示验证。
- 56 个扫描帧的 10 个来源家族全部保留。新增排重底账共有 383 个记录：273 个原案例整体及 110 个可能的组分作用域；后者是待审假设，不是新增正式案例。所有化学作用域仍需审阅，同组成的未决原库记录会阻止自动新增计数。纯粹旋转、基组、自旋/电荷表示及环境组分变更不能自动获得新增分子信用。
- OLD-099 的原输入标题说明固定几何，而质量参考登记为平衡候选，作为用途差异保留到 REF-001 全批结束后统一审阅；不在当前批次修改作业或据此判为输入错误。
- 正式完整通过仍为 0，新增外部分子仍为 0，生产科学问题未关闭。可读证据：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\input-structure-identity\GI-001\post-batch-review-v1\原库结构取证与去重底账.md`；距离复核：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\input-structure-identity\GI-001\distance-control-followup-v1\原库结构距离判据复核.md`；机读结果：`F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\input-structure-identity\GI-001\post-batch-review-v1\review.json`。


## 2026-09-06 拓扑通解实施前约束补充

- 已建立 general-fix-contracts/topology-v1，记录 8 份源码/既有测试身份和 6 个实施单元，状态 prepared_not_applied。
- TOPO-001 分别处理上游方向快捷判据、下游多通道强制 Spiro 和紧凑集合消费者；保留真实结构见证、方向通道与电子证据的不同作用域。
- IUPAC 的一般螺连接与自由螺连接不同：割点/循环块只作为自由连接结构证据，不能凭没有割点排除所有带附加桥联系的螺结构。GI-001 环块参考范围据此明确，原始取证和控制结果不改写。
- TOPO-002 保持最新裁决：OLD-013 主壳 4Cl 加另一 Re、CN5；修订远侧联系的次级显示层，不恢复已撤回的第一壳/CN 算错结论。
- 本补充没有关闭科学问题、改变正在运行的 REF-001，也没有新增正式通过或外部分子。详细契约：F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\general-fix-contracts\topology-v1\拓扑与次级联系通解实施契约.md。


## 2026-09-06 REF 暂停交接与首批通解实施检查点

用户最新指示优先于旧契约的 REF 前置门槛：先推进全部不依赖 REF-001 的修复。
REF-001 已停止并核对所属进程及运行锁释放，状态为 13 个候选采集完成、
2 个中断保留（OLD-013/015）、258 个未开始。原始阶段记录、检查点和中间文件
保留；当前暂停状态由单独覆盖记录说明。没有将候选完成计为物理验收通过。

GitHub 交接已完成：
[资料包](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906)，
[固定源码和 273 组实际输入](https://github.com/O1dDing/CUDA-Orbital-Visualisation/tree/51c9815b28f26e796be7892d552ccf8ad86b8785/validation/paused-20260906)。
23 个分卷共 9,837,529,264 字节，保留 78,481 个文件。服务器大小/摘要、公开
清单回读、远端源码对象和本地全量恢复校验均通过；恢复检查验证 273 组输入、
28 个已完成阶段，未启动新 Gaussian 计算。

生产数值修复源码 `a9bf30a18a16c16e876abec133691b911696da3f` 已构建普通版和
验证版：统一 Gaussian AO 相位/排列，CPU 参数和 MO 保留 double，GPU 显式
上传；S 从 s–g 基组解析积分得到，支持矩形 MO，全自旋 Gram/谱/秩/残差有
明确状态；producer S 独立保留核验；CHEM-001 累加器从零开始。
NUM-001/002/003、CHEM-001 状态为“已实施并通过组件检查，待全库验证关闭”。
合理的价层和可读性政策尚未据此全部验收，不误关其他问题。

ON/OFF 各 33 项组件测试通过，0 失败、0 跳过。首个原库案例 OLD-001 的
15 个 MO 采样和 4 个关键轨道完整网格通过既定独立门槛，27 张截图、197 条
实际控件操作及四种导出完成，0 操作失败。包含普通版取证和整理约 10.330 秒，
独立检查主体约 1.951 秒；这不是包含 Gaussian 补算和完整科学/UI 裁决的总耗时。
compact 截图未完整露出紧凑图，保留图像覆盖观察，未给出全 UI 通过结论。
最初缺少 Pillow 的采集环境失败保留在 NUM-FIX-001-pilot，未计作通过。

当前全 273 轮 NUM-FIX-003-original 使用固定程序和输入，四进程、每进程
3 核/24 GiB，整个任务受 12 核/128 GiB 限制。全部案例终止后统一比较和分析，
轮内不改程序、输入或门槛。正式完整通过 0，新增正式外部分子 0。
最新机读进度位于 outputs/general-fixes-20260906/implementation-progress.json。


## NUM-FIX-003 原库 273 全量取证完成


同一程序、输入和门槛下，273 个原案例全部取证结束：
{'complete': 273, 'completed_with_gaps': 0, 'error': 0}。共 31,871 个 MO
实际纹理、12,191 张截图、
1,214 份前线完整纹理。
取证墙钟 2100.375 秒（约 35 分钟）。
这些是取证完整性结果，独立数值、实际图像及科学解释继续复核；正式完整通过仍为 0。
固定轮次身份 `4c6a0986839ecaf985be304d66b75e89dec54288400c1839423a3c681348aa45`。错误/缺失没有被移出 273 清单，
REF-001 继续暂停。本次已有 g 壳层控制输入复查不计入新增外部分子。


## NUM-FIX-003 数值复核完成与下一版对称性修复


REF-001 按用户指示保持暂停：13 个候选采集完成，2 个中断保留，258 个未开始。
[完整 REF 恢复交接](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906)
已上传并验证。该资料包不是最终软件发布。

生产源码 `a9bf30a18a16c16e876abec133691b911696da3f` 的 NUM-FIX-003 原库轮已完成：

- 273/273 完整取证，273/273 通过本轮独立数值子集检查，0 失败、0 缺证。
- 31,871 个 MO 均有 8,179 个唯一网格采样点；1,214 份前线完整纹理通过固定门槛。
- 普通/验证构建记录的科学数据一致；真实验证 EXE 的输入诊断与同源读取结果一致。
- 独立 AO 积分最大绝对差 1.66e-14；最大 MO Gram 残差 1.03e-6。
- 全 MO 采样最大 NRMS 7.22e-5（固定门槛 1e-4）；完整前线纹理最大 NRMS 1.10e-5。
- 400 AO / 396 MO 的矩形案例也通过独立基组积分、Gram 与电子数检查。
- 另有纯球谐 220 MO、笛卡尔 280 MO 的 g 壳层控制通过。它们不计作外部新分子。

原库取证墙钟 2,100.375 秒，独立复核 1,086.711 秒。单案例取证中位数 24.158 秒、
95 分位 58.726 秒；独立复核中位数 12.606 秒、95 分位 33.903 秒。
原始萘 OLD-015 分别为 38.862 秒和 23.867 秒。未包含暂停的 Gaussian 参考补算与
尚未完成的图像、科学解释裁决。网格比较使用已记录的 CUDA 浮点坐标插值模型；
采样和有限完整网格通过不等于连续全空间证明。

NUM-001/002/003、CHEM-001 已实施并得到上述范围的验证支持。完整产品问题台账继续
保留其验证版本与适用边界。正式完整通过案例仍为 **0**，正式外部分子仍为 **0**。
全部存储密度字段、所有规定不变量和剩余 UI/化学问题没有被算入已经通过的范围。

后续源码 `94447b057be12f3126f6cbef8743bb894c622144` 已修正偶数阶 Dnh 的 g/u 与
B 类标签、Gaussian 线性缩写排版及空点群字段读取，并分开 detected/used 来源。
保留原始点群行、缺项/无效状态、Link1/行号，以及派生 MO 子空间、坐标轴约定和残差。
ON/OFF 各 33 项组件检查通过，0 失败、0 跳过；这一新版本的完整原库回归仍待执行。

Dnh 分类依据对称操作的实际字符及反演关系，B1/B2 绑定所记录的 C2 轴约定。
[原作者 D6h 字符表及轴约定](https://gernot-katzers-spice-pages.com/character_tables/D6h.html)
用于交叉核对；线性缩写与 [cclib Gaussian 原始解析器](https://github.com/cclib/cclib/blob/master/cclib/parser/gaussianparser.py)
交叉核对。已有 producer 标签原文保持独立，符号正规化只作用于排版。

此目录保存本轮清单、逐案例数值结论、误差极值、耗时、参考环境身份、脚本及组件结果。
大型新取证纹理/截图/独立网格保存在 status.json 所列本地目录；本次紧凑进度提交没有
宣称它们全部上传。REF 交接资料的完整上传状态不受影响。

继续工作顺序为：剩余全局/局部及 MO/显示组作用域 → 对称推导与混合子空间 → 环拓扑
见证和联系层次 → 配体先验及能隙类型 → 界面/导出同一快照 → 新版本全 273 复验。
REF 与外部 Gaussian 计算保持暂停，等待用户以后分批恢复。


## NUM-FIX-004 当前视图快照与全量回归

源码构建身份 `70e72ef9d4c956b374a36efb4a478489127bca90`，包含上一阶段对称性修复及当前视图快照通解。验证版与普通版均已编译，34 项组件检查各自全部通过。

- EXP-001 已实施：真实导出按钮传递已经绘制的不可变快照；原始 MO、决定行选择的 anchor 和当前实际检查的 MO 分开记录。后续控件变化不会回写旧快照。
- PNG/SVG 使用完整实际成员及电子箭头，β 对应关系和选中成员保留；JSON/CSV 记录实际成员原始值、行布局值、成员能量范围及可见行映射。机读浮点输出保留 double 往返精度。
- COV_VALIDATION 1 保持兼容；新增只控制证据文件名的 export-name 指令，导出仍经过真实按钮。当前轮分别保留至少 7 个视图状态，含 β 对应成员的案例额外取证。
- 首个 OLD-001 完整取证 11.7034204 秒，独立数值检查 1.9584581 秒；7 组导出视图检查全部通过。该项仅通过相应子集，不等于完整化学或界面验收。
- 已查看真实 PNG 和 compact-top 截图。能级图仍不能在当前卡片内一次完整显示，摘要也有裁切；保留 UI 的相关裁决。绘制调用、选中 PNG 线条探针和快照身份检查不能代替全部文字/像素验收。
- 已启动原 273 全量取证，轮身份 `79a352a6ef9de3413505f2c8c25c06561d6e5da2086240a8f836cd4da3b30693`。所有案例终止后再做本轮统一复核；批次中不修改生产代码、输入或阈值。
- REF-001 继续暂停；正式完整通过 0，新增外部分子 0。既有 NUM-FIX-003 的 273 数值通过证据仍属于 a9bf30a，不能归到本版。

机读进度：`F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906\view-snapshot-progress.json`。活动轮进度：`F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906\NUM-FIX-004-original\progress.json`。


## NUM-FIX-004 全 273 复核完成与重叠布局问题

程序构建 `70e72ef9d4c956b374a36efb4a478489127bca90`，轮身份 `79a352a6ef9de3413505f2c8c25c06561d6e5da2086240a8f836cd4da3b30693`。全部取证及复核进程已结束；无 Gaussian 新作业。

- 273/273 取证完整，无崩溃、超时或证据缺项。273/273 通过本轮独立数值子集，31,871 个 MO 的实际纹理采样及 1,214 份完整前线纹理均通过。普通版/验证版的完整科学输出与实际 EXE 数值证据在全库相符。
- 多状态视图复核完成 273/273。1,943 组四格式导出均通过快照身份、实际成员、原始能量、电子箭头和导出帧关联检查；其中 4 组另有选中 PNG 内部线条像素失败。因此视图案例状态为 269 通过、4 未通过，不能写成全库导出已全部通过。
- 新增 **EXP-002：近邻能级的线条及点击区域重叠**。OLD-071、OLD-175、OLD-176、OLD-249 的 export-linear 状态中，后绘制能级覆盖了先绘制的选中线条内部。OLD-071 的蓝色外缘仍可见，不将它描述成整个选中标记消失。其四行图仅有两处可分辨的线条位置，近邻行的命中区域也相互重叠。
- 根因：屏幕和导出只把能量四舍五入到 1e-12 Ha 后完全相同的行横向错开。实际相差 1e-5 Ha 左右的行在线性轴上可能远小于一个像素；它们既不会归入同桶，也不能在原横向位置分辨。新修复将根据实际绘制/命中范围分配横向通道，保持原始能量、行分组和纵坐标，不把近邻强制当作严格简并。该布局修复尚未实施/验证。
- EXP-001 的“导出重新构建而丢失当前视图”根因在本版本的全库多状态检查中未再出现；保留 EXP-002 和完整图像/文字验收，按子项裁决。SYM-004/005、POL、TOPO、其余 UI 与物理参考问题仍保留。
- 总取证墙钟 2182.597 秒；独立数值复核 1090.602 秒；视图复核 24.977 秒。单案例取证中位数 25.455 秒、P95 60.295 秒；数值复核中位数 12.481 秒、P95 34.684 秒。
- OLD-015 萘本轮取证 41.415758 秒，独立数值复核 24.582416 秒，合计 65.998174 秒。这不包含 Gaussian 参考计算及尚未完成的化学/完整界面验收，也没有把批量视图复核的 24.977 秒均分成未经测量的单案例耗时。
- REF-001 按用户指示暂停；正式完整通过 0，新增外部分子 0。所有失败与冻结证据保留。对话服务曾出现的 Bad Request 不属于本轮 COV/Gaussian 科学失败。

证据：`original-numeric-summary`、`original-numeric-case-reports`、`original-view-review`、`original-unified-review`。其中每个视图报告均记录消耗证据的哈希。仍需完成后续不依赖 REF 的通解修复及全库循环。


## NUM-FIX-005 近邻能级布局通解与原 273 新一轮

构建源码 `e1e118d6dc34839932c2f28a0fd8379a40e4f3b4`。上一轮完整分析保留的 EXP-002 已实施共同根因修复：按实际绘制与点击范围进行横向通道分配，屏幕与 PNG/SVG 使用同一布局算法。能量纵坐标、原始成员和真实能量均未改变；不会把像素接近当作严格简并。需要额外宽度时屏幕画布支持真实横向滚动，导出画布扩大，PNG/SVG 尺寸一致。

普通版与验证版各 35 项组件检查通过，无跳过。新增检查从实际 ImGui 绘制顶点定位相邻 MO，分别通过鼠标点击验证选择，并用真实横向滚轮检查溢出通道。保留完整文字/像素验收边界。

首个 OLD-001 试运行已完成：单案例取证 12.569970 秒，独立数值复核 1.992701 秒；7 组多格式导出通过，全部线性轴成员点击通过。原库 NUM-FIX-005 已冻结并启动，全部 273 终止后自动进行整库数值与视图复核。本轮增加逐成员线性点击、实际命中范围重叠、导出线条重叠与画布尺寸检查，不降低原数值或像素阈值。

整库结果未出，EXP-002 尚未关闭。SYM-004/005、局部 S 度量投影、TOPO、POL、其余 UI 和独立物理参考仍待处理。REF-001 保持暂停，正式完整通过 0，外部分子 0。

活动目录：`F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906\NUM-FIX-005-original`。`progress.json` 保存取证进度；`pipeline-progress.json` 保存顺序复核状态。不得仅因对话服务中断而重复启动已有进程。此目录只包含紧凑恢复记录；不声称本轮全部纹理/截图已上传。


## NUM-FIX-005 完整终态裁决与 NAV-001

冻结生产版本 e1e118d6dc34839932c2f28a0fd8379a40e4f3b4。全 273 取证于 2026-09-06 23:57:59 结束，249 完整、24 有交互缺项、0 运行错误；统一数值及视图复核于 2026-09-07 00:17:05 结束。随后对话的额度中断没有中断后台复核，本次恢复不重复运行这轮。

- 数值子项：273 全通过，31,871 个实际渲染纹理的同点样本、1,214 个完整前线网格；门槛未改变。
- 多格式导出：1,943 组全部通过，含来源/快照/成员/能量/箭头、PNG 选中线条、PNG/SVG 尺寸、线条及实际命中范围无重叠。
- 逐成员线性轴点击：2,699 通过、72 失败；对应 24 案例视图子项失败，其余 249 通过。72 个目标的点击与悬停合计 144 次失败，缺文件和完整性错误均为 0。
- EXP-002 的近邻线条遮挡及命中区域重叠子项在本轮检查范围内已修复。不能据此宣称所有图像文字、交互或化学结论已通过。
- 新增 NAV-001：嵌套滚动视图中的横向目标不可达。全部 72 个失败目标的纵坐标已在实际绘制裁切范围内，但目标不可见，实际选择未切换。源码 seek 会在横向裁切时仍先对齐纵向父子窗口；横向滚动条还可能被外层卡片裁掉。同帧鼠标移动/滚轮与逐帧 ClearEventsQueue 的组合也需按 ImGui 的事件分帧机制修正。下一版保留真实鼠标及滚轮路径，不直接设置滚动值或选中 MO。

计时：取证批次 2402.582898 秒、独立数值复核 1095.476959 秒、视图复核 27.611610 秒。取证中位数 28.243681 秒、95 分位 72.770153 秒；OLD-015 萘取证 52.896140 秒、数值复核 24.964995 秒，合计 77.861135 秒，不含 Gaussian 或完整化学验收。

2026-09-08 恢复时发现用户另在 F:/CalChem/COV/Resume 启动了独立增量参考流程：workers=2，单作业 4 核/24 GB，OLD-018 initial 仍活动，PAUSE.json 请求阶段边界暂停。其旧 progress.json 是恢复拷贝，不代表当前进程。主任务不终止、重启或修改该流程，保留 8 核/48 GB 的预算，本任务最多使用余下 4 核/64 GB。冻结历史报告中的 REF 暂停描述只代表历史批次时点，不能当作当前无 Gaussian 的证据。

正式完整通过 0、外部正式分子 0。SYM、局部投影、TOPO、POL、剩余 UI、物理参考等未完成项目继续保留。当前保存的是完整紧凑报告及其校验值；本轮大体积原始纹理和全部截图未上传。


## NUM-FIX-006 导航通解已实施，开始全 273 复验

源码版本 `21036b2990ef86dcdc419c65cdebd2a03bf4a12b`，全库冻结身份 `389d5ae5700ec4ab04b38b814b533b47b7cf1c646e0a0d09dd7829d2b4d7403b`。验证版和普通版各 36 项组件检查全部通过。鼠标坐标按 ImGui 的整数像素处理，移动和滚轮分帧，目标按实际父子窗口可见范围逐轴导航；没有直接改滚动值或选中轨道。

指定首例 OLD-001（H₃⁺）：取证 17.177897 秒、独立数值单案例 2.054437 秒、视图检查 0.067088 秒；15 个实际 MO 纹理样本、4 个完整前线网格和 7 组导出快照均通过。这些时间不含 Gaussian、独立环境身份核验开销和全部物理验收。

额外诊断 OLD-013：216 个 MO 纹理样本、4 个完整前线网格、7 组导出快照及全部 21 次线性轴成员点击通过。上一轮失败的零基成员 [57, 59, 60, 61] 均通过真实控件选中；所附 MO 58 截图对应内部成员 57，可见实际选中线和提示文字。取证 50.757531 秒，独立数值单案例 37.396772 秒。仍可见长提示框底部裁切，UI-004 继续保留，不能将导航通过扩大为全部界面通过。

全 273 采集和其后的自动复核已启动，只有采集监督进程退出并确认全部案例终止，流水线才进行独立数值与视图复核。全库轮内固定程序、输入、判据和取证器；不提前关闭 NAV-001。当前机器进度另存原子更新文件，不将部分采集状态当作全轮结论。

本任务使用 2 个进程、每个 2 核/24 GiB，并由共同父进程限制在 4 核/64 GiB 内。为用户独立恢复的 Gaussian 队列继续预留 8 核/48 GiB，不改变其阶段边界暂停请求。旧冻结报告中的 REF 暂停字段是历史固定文案，不代表现在没有 Gaussian 进程；本次运行状态以此处明确的独立队列管理边界为准。

正式完整通过仍为 0，正式外部分子仍为 0。SYM、局部度量投影、POL、TOPO、密度补证、其余 UI、物理参考和外部正式集合均继续推进。


## NUM-FIX-006 全 273 终态裁决：NAV-001 已修复

固定生产版本 `21036b2990ef86dcdc419c65cdebd2a03bf4a12b` 的 273 个案例均完整取证，零崩溃/超时/证据缺失；独立数值子集 273/273，31,871 个全 MO 同点纹理及 1,214 个完整前线网格均满足本轮冻结门槛。1,943 组视图/导出快照检查全部通过。

逐条核对全部 2,771 个线性轴成员点击及上一轮相同目标集合。上一轮 24 个案例中的全部 72 个失败目标现均通过，目标覆盖没有减少。NAV-001 在这套实际鼠标/滚轮、成员选择与对应帧证据范围内记为已修复。没有通过直接赋值滚动量或所选 MO 来代替真实控件。

UI-004 长提示框底部裁切继续保留；视图子集通过不等于所有字形、布局或全部像素已正确。密度来源、作用域、拓扑、化学解释、全不变量及物理参考仍需后续完整闭环，正式案例完整通过数仍为 0。用户管理的 Gaussian 作业状态独立保留，不沿用旧摘要中的“全部暂停”作为当前控制状态。

本轮采集墙钟 9,727.501 秒，独立数值复核 2,387.533 秒，视图复核 31.461 秒。案例采集中位 79.584 秒、P95 134.525 秒；这些耗时不包含 Gaussian 补算，也不是完整物理验收耗时。全量逐案例与 OLD-015 的本轮计时均保存在统一 JSON 摘要。


## 密度旧生产库 20 项完整对照：新增 DEN-001 至 DEN-005

生产版本 `21036b2990ef86dcdc419c65cdebd2a03bf4a12b` 的 20 项输入全部进入终态后统一审查：8 项满足固定要求，12 项失败或缺证。共享整数占据在 MO 重排后的 Q 最大变化为 1.2614725579091308，远超事先固定的矩阵界限 2e-12。原有正确的显式两自旋、已知零 Q、虚轨道截断、显式分数占据和系数不完整拒绝等 8 项行为列为保留正例。

- **DEN-001（新增，数据有效性与可用状态）**：缺失或非法 Occup 仍产生带 Derived 来源的密度，零值、缺项和失败状态混淆。 通解：保留逐 MO 字段来源，严格检查普通占据模型的范围和有限性；P/Q 分别保存矩阵、状态和原因。 对照：DEN-C01, DEN-C11, DEN-C12, DEN-C13, DEN-C14。
- **DEN-002（新增，轨道身份与占据解释）**：共享整数占据的未配对方向按列表索引重新填充；同一 MO 集合重排使 Q 改变。 通解：规范 FCHK 占据绑定不可变来源身份；Molden 保持逐 MO Occup 并明确共享整数单行列式假设。 对照：DEN-C02, DEN-C03。
- **DEN-003（新增，自旋来源与模型边界）**：共享分数占据或缺失 Spin 时，未知 Q 被表示为可用矩阵或带 Derived 的空矩阵。 通解：显式 α/β 与共享空间分别处理；共享分数占据不自动推断自旋分区，已知 P 保留。 对照：DEN-C04, DEN-C05, DEN-C10。
- **DEN-004（新增，占据子空间完整性）**：缺少占据方向的原始 MO 重建返回完整外观的 P/Q；C19 的输入矩阵正确但原始重建仍错误。 通解：用独立输入计数和来源身份检查占据覆盖；实际 producer P/Q 与原始 MO 重建分开保存。 对照：DEN-C18, DEN-C19。
- **DEN-005（新增，派生分析的输入依赖）**：P 或 Q 不可用或不正确时，完整 Mayer 结果仍被标记 Derived。 通解：完整 Mayer 分析要求两个有效密度；缺 Q 时不省略自旋项后沿用相同名称。 对照：DEN-C01, DEN-C02, DEN-C03, DEN-C04, DEN-C05, DEN-C10, DEN-C11, DEN-C12, DEN-C13, DEN-C14, DEN-C18。

这些问题允许重复关联同一输入，不把 12 个失败控制虚增为新的分子数。修复源码 `e2e835a70229b94501127ad336e0ddbc3ea69df3` 已实施，正在编译及组件检查；尚未完成同一检查器的修复复验，以上问题此时均不关闭。完整病例通过数为 0。基线矩阵、原始输出、进程退出状态及检查器字节全部保留。


## 密度通解固定对照通过、详情试运行缺口与 50 个外部结构去重

密度通解 e2e835a7 已用同一冻结检查器复验普通版 20 项、验证版 20 项。全部 40 项先进入终态，再统一比较；40/40 满足原固定要求，20 组 ON/OFF 实际输出完全一致，原矩阵界限 2e-12 未改变。DEN-001 至 DEN-005 更新为“已实施通解且固定对照通过，待完整 273 回归关闭”，没有把数学控制算作新增分子。

UI-004 已实施短提示与可滚动详情窗口。NUM-FIX-007-pilot 的 OLD-001 全部采集结束，取证 24.196287600003416 秒；数值及导出子集通过，但原生操作有 2 条失败，整次试运行保留为有缺口。真实关闭已生效，自动操作器却继续要求已消失的关闭按钮存在，随后误报目标未绘制。新增 EVID-ACT-001：真实关闭动作的完成判定不应依赖目标继续存在；20a3263 已修复，待后续实际试运行及全库检查。

新增 EVID-CLIP-001：可点击的部分可见目标不能证明说明文字完整可读。已查看实际详情顶部和底部 PNG；第一版末尾说明仍部分裁切。后续计划先通过真实滚轮到达底部关闭按钮再截图，并保存真实 ImGui 裁切矩形；独立检查器核对完整矩形及全部失败动作。该问题保留原反例，不用自报 visible=true 关闭裁切问题。UI-004 此时仍不关闭。

8b25227 增加普通/验证构建共同的实际 P/Q 序列化及独立矩阵检查：从同一 FCHK 的 IOData 系数/占据、原始 producer SCF 密度和独立 AO 约定建立参考；输入矩阵变换要求精确一致，MO 重建只允许由运算次数决定的双精度误差，不借用十进制打印精度。新检查器的 10 项正反例通过，旧试运行实测密度的 17 项诊断通过；这是冻结下一轮前的检查器准备，不改写旧轮判据。

50 个外部起始结构的全部同分子式冲突已解决：噁唑/异噁唑的 O–N 连接数不同；苄基自由基的六元环加支链与 OLD-097 七元环的碳连接度序列不同。各自已冻结的四档距离检查与作者/PubChem 连接证据一致。全部 50 个彼此不同且与原库不同；尚未完成参考质量和电子态验收，正式外部案例仍为 0。

Gaussian 继续由用户的 Resume 流程独立管理。本检查点没有启动、重启、终止该队列，也没有修改其暂停文件。COV 工作沿用物理核槽 8–11，最多 4 核和 64 GiB 的进程树限制。


## NUM-FIX-007C：试运行与长提示诊断通过，已冻结启动全 273

源码 `aaabc3edb13c2768f6bc9063ee933a8881e7f57f`，取证/复核源码 `93e294b0bca5d920494cfbb1ab76f4c218aaf898`。ON/OFF 各 37 项组件测试通过，0 失败、0 跳过。NUM-FIX-007 和 007B 的试运行缺口继续保留；修复关闭动作完成判定并让悬停通过真实滚轮完整显示可容纳的目标后，007C 的首个案例及 OLD-013 长提示诊断全部已冻结子检查通过。

- OLD-001：取证 26.637406 秒，独立数值检查本体 1.832030 秒，视图复核 0.094225 秒；15 个 MO、4 个完整前线网格、7 组导出、4 张详情控件实帧通过。
- OLD-013：取证 50.162333 秒，独立数值检查本体 30.169827 秒，视图复核 0.221731 秒；216 个 MO、4 个完整前线网格、7 组导出、4 张详情控件实帧通过。

全库轮身份 `475be3b2459c105da3686e3065205e4e7bcfde895cacdb6d0457f036f214b1e3`。273 个输入及其旧编号全部冻结，2 个 COV 进程并行，每个 2 线程；整体 4 个物理核槽、64 GiB 进程树上限。只有全体采集终止且监督进程退出后才开始冻结的独立数值、密度和视图复核。不会重跑 NUM-FIX-005，也不会在这轮中改变算法、判据或输入。

实际查看了首个案例的详情底部/重新打开帧以及 OLD-013 的 MO 58 短提示帧；指定内容完整可见。此范围不代表全部字形、科学解释或所有像素已通过。DEN-001 至 DEN-005、UI-004 及本次取证问题均等待全库回归裁决；全部正式案例完整通过数仍为 0。Gaussian 由用户独立队列管理，未执行控制操作。
