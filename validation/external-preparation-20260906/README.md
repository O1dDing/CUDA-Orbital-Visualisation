# 外部集合准备：67 个候选，尚未冻结正式集合

54 个候选的公开身份字段经电荷和元素计数核对后相符，其中 51 个的元素计数与原库全部 273 案例不同。另有 9 个名称查询无匹配、1 个数据库字段与目标不符、3 个按原始结构论文继续核对。全部原始响应、错误响应和复核过程保留。正式外部分子数仍为 **0**；没有启动 Gaussian，也没有形成可验收的外部 FCHK。

查询按 [PubChem PUG REST 官方接口](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)顺序进行，并限制请求速率。保存 CID、InChIKey、连接 SMILES、分子式、电荷、响应原文和 SHA-256。第一版读取器未支持 PubChem 分子式尾部的电荷表示，第二次复核直接重读原始响应并显式比较电荷，未改写外部数据。

同分子式不等于同一结构。cubane/COT、pentatetraene/spiropentadiene、benzyl/tropylium 分别保留原库碰撞记录，正式入库前还须核对图结构。54 个相符候选的连接 InChIKey 首段之间没有重复；金属配合物的数据库标准化可能拆开配体，此结果不能证明配位连接、几何或电子态。

六甲基钨的名称响应给出了断开的甲基组分及 −6 电荷，与目标中性 W(CH3)6 不符，已作为输入反例保留，不能直接用于 Gaussian。Ni(PMe3)4 名称查询同时返回质子化与中性记录，逐条核对后仅保留中性匹配的身份候选。[Ni(PF3)4 气相结构原始研究](https://pubs.rsc.org/en/content/articlelanding/1970/c2/c29700000595)和 [Rh/PMe3 配合物原始结构研究](https://pubs.rsc.org/en/content/articlelanding/1980/dt/dt9800000511)提供后续坐标核查入口，当前尚未取得并验收这些坐标。

甲基过氧自由基采用 [NIST CCCBDB 的身份标识](https://cccbdb.nist.gov/exp2x.asp?casno=2143580&charge=0)补查；硝普盐阴离子采用 [PubChem 对应阴离子条目](https://pubchem.ncbi.nlm.nih.gov/compound/11963622)补查。保留最初的查无结果或中性错误匹配，不能让名称别名静默改变目标电荷。

正式集合还需补充弱复合物/环境模型，例如氨—氟化氢、HCN—水、甲酸二聚体、甲醇—水及 CO2—水；这五项目前是用途候选，尚未另计入 67 项身份清单。新增分子中的纯/笛卡尔 g 表示作为同一分子的输入变体覆盖，不额外增加分子计数。

后续顺序：核对真实结构或适用初始几何 → 明确电荷、电子态及环境用途 → 对原库和新增集合进行化学图去重 → 明确冻结成员、方法和判据 → 在用户解除 Gaussian 暂停并满足前置条件后计算及正式验证。准备中的失败条目继续保留；不得把身份匹配数当作计算完成或 COV 通过数。

| 准备编号 | 候选 | 分子式 / 电荷 | 当前身份状态 | 用途 |
|---|---|---|---|---|
| PREP-001 | silane | SiH4 / 0 | 身份字段相符；几何待建 | 四面体及 Si 价层 |
| PREP-002 | germane | GeH4 / 0 | 身份字段相符；几何待建 | 重主族及 ECP 对照 |
| PREP-003 | phosphine | PH3 / 0 | 身份字段相符；几何待建 | P 孤对与后续配位对照 |
| PREP-004 | phosphorus trifluoride | PF3 / 0 | 身份字段相符；几何待建 | P 配体先验与轨道证据 |
| PREP-005 | phosphorus trichloride | PCl3 / 0 | 身份字段相符；几何待建 | 不同卤素与 P 价层 |
| PREP-006 | sulfur tetrafluoride | SF4 / 0 | 身份字段相符；几何待建 | 畸变配位及多中心 |
| PREP-007 | chlorine trifluoride | ClF3 / 0 | 身份字段相符；几何待建 | T 形与孤对 |
| PREP-008 | bromine pentafluoride | BrF5 / 0 | 身份字段相符；几何待建 | 方锥与重元素 |
| PREP-009 | iodine heptafluoride | IF7 / 0 | 身份字段相符；几何待建 | 七配位与重元素 |
| PREP-010 | oxygen difluoride | OF2 / 0 | 身份字段相符；几何待建 | 弯曲分子及键极性 |
| PREP-011 | sulfur trioxide | SO3 / 0 | 身份字段相符；几何待建 | 平面与离域 |
| PREP-012 | carbon disulfide | CS2 / 0 | 身份字段相符；几何待建 | 线性群及多个 π 通道 |
| PREP-013 | carbonyl sulfide | COS / 0 | 身份字段相符；几何待建 | 异核线性及局部不等价 |
| PREP-014 | nitrous oxide | N2O / 0 | 身份字段相符；几何待建 | 异核线性与多中心 |
| PREP-015 | hydrazine | N2H4 / 0 | 身份字段相符；几何待建 | 非平面与孤对耦合 |
| PREP-016 | nitrosyl chloride | ClNO / 0 | 身份字段相符；几何待建 | 混合键极性 |
| PREP-017 | thionyl chloride | Cl2OS / 0 | 身份字段相符；几何待建 | S/O/Cl 价层与局部几何 |
| PREP-018 | borazine | B3H6N3 / 0 | 身份字段相符；几何待建 | 杂原子环及离域 |
| PREP-019 | phosphinine | C5H5P / 0 | 身份字段相符；几何待建 | 环内 P 与多方向轨道 |
| PREP-020 | pyrazine | C4H4N2 / 0 | 身份字段相符；几何待建 | 多个 N 与全局对称性 |
| PREP-021 | 1,3,5-triazine | C3H3N3 / 0 | 身份字段相符；几何待建 | 三重 N 与简并空间 |
| PREP-022 | oxazole | C3H3NO / 0 | 身份字段相符；几何待建 | N/O 混合 |
| PREP-023 | isoxazole | C3H3NO / 0 | 身份字段相符；几何待建 | 同分子式异构图去重 |
| PREP-024 | thiazole | C3H3NS / 0 | 身份字段相符；几何待建 | N/S 混合 |
| PREP-025 | indole | C8H7N / 0 | 身份字段相符；几何待建 | 稠环与孤对 |
| PREP-026 | benzofuran | C8H6O / 0 | 身份字段相符；几何待建 | 稠合 O 杂环 |
| PREP-027 | benzothiophene | C8H6S / 0 | 身份字段相符；几何待建 | 稠合 S 杂环 |
| PREP-028 | pyrene | C16H10 / 0 | 身份字段相符；几何待建 | 多环 π 网络 |
| PREP-029 | perylene | C20H12 / 0 | 身份字段相符；几何待建 | 较大 π 空间 |
| PREP-030 | cubane | C8H8 / 0 | 身份字段相符；几何待建 | 高对称 σ 骨架及 COT 图去重 |
| PREP-031 | spiro[3.3]heptane | C7H12 / 0 | 身份字段相符；几何待建 | 真实螺连接正例 |
| PREP-032 | bicyclo[1.1.1]pentane | C5H8 / 0 | 身份字段相符；几何待建 | 桥联骨架与非螺反例 |
| PREP-033 | norbornadiene | C7H8 / 0 | 身份字段相符；几何待建 | 空间分离双键及桥骨架 |
| PREP-034 | pentatetraene | C5H4 / 0 | 身份字段相符；几何待建 | 多通道及非螺反例 |
| PREP-035 | malonaldehyde | C3H4O2 / 0 | 身份字段相符；几何待建 | 环境及构象需另行确定 |
| PREP-036 | acetylacetone | C5H8O2 / 0 | 身份字段相符；几何待建 | 分子内氢键与互变异构 |
| PREP-037 | biphenyl | C12H10 / 0 | 身份字段相符；几何待建 | 非共面芳环 |
| PREP-038 | 1,4-dioxane | C4H8O2 / 0 | 身份字段相符；几何待建 | 孤对与构象 |
| PREP-039 | acetic acid | C2H4O2 / 0 | 身份字段相符；几何待建 | 单体及后续弱复合物对照 |
| PREP-040 | dimethyl sulfoxide | C2H6OS / 0 | 身份字段相符；几何待建 | S/O 极性与后续配位 |
| PREP-041 | TEMPO | C9H18NO / 0 | 身份字段相符；几何待建 | 开放壳层及 α/β |
| PREP-042 | phenoxy radical | C6H5O / 0 | 身份字段相符；几何待建 | 离域自由基 |
| PREP-043 | benzyl radical | C7H7 / 0 | 身份字段相符；几何待建 | 与环状 C7H7 图去重 |
| PREP-044 | ethynyl radical | C2H / 0 | 身份字段相符；几何待建 | 线性开放壳层 |
| PREP-045 | methylperoxy radical | CH3O2 / 0 | 身份字段相符；几何待建 | 过氧与自旋 |
| PREP-046 | hydroperoxyl radical | HO2 / 0 | 身份字段相符；几何待建 | 小型开放壳层 |
| PREP-047 | tetrakis(trifluorophosphine)nickel | F12NiP4 / 0 | 身份字段相符；几何待建 | P 先验与 π 受体轨道 |
| PREP-048 | tetrakis(trimethylphosphine)nickel | C12H36NiP4 / 0 | 身份字段相符；几何待建 | 候选名称及真实结构待核对 |
| PREP-049 | pentacarbonyl(trimethylphosphine)chromium | C8H9CrO5P / 0 | 名称查询无匹配，保留 | CO/PMe3 |
| PREP-050 | pentacarbonyl(trifluorophosphine)chromium | C5CrF3O5P / 0 | 名称查询无匹配，保留 | CO/PF3 |
| PREP-051 | cisplatin | Cl2H6N2Pt / 0 | 身份字段相符；几何待建 | 顺式 Pt 配合物及局部/整体 |
| PREP-052 | dichlorobis(trimethylphosphine)platinum(II) | C6H18Cl2P2Pt / 0 | 名称查询无匹配，保留 | 顺反结构不得凭名称选定 |
| PREP-053 | chlorotris(trimethylphosphine)rhodium(I) | C9H27ClP3Rh / 0 | 名称查询无匹配，保留 | 实验畸变及局部群 |
| PREP-054 | carbonylchlorobis(trimethylphosphine)rhodium(I) | C7H18ClOP2Rh / 0 | 名称查询无匹配，保留 | CO/PMe3/Cl |
| PREP-055 | bis(ethylenediamine)copper(II) | C4H16CuN4 / 2 | 身份字段相符；几何待建 | 有限环境与电子态待建 |
| PREP-056 | bis(ethylenediamine)zinc(II) | C4H16N4Zn / 2 | 名称查询无匹配，保留 | 螯合几何 |
| PREP-057 | tris(ethylenediamine)cobalt(III) | C6H24CoN6 / 3 | 身份字段相符；几何待建 | 手性及不同作用域 |
| PREP-058 | nitroprusside ion | C5FeN6O / -2 | 身份字段相符；几何待建 | NO/CN 非无辜配体 |
| PREP-059 | permanganate | MnO4 / -1 | 身份字段相符；几何待建 | 高氧化态与简并 |
| PREP-060 | molybdate | MoO4 / -2 | 身份字段相符；几何待建 | 与 MnO4 不同电子结构 |
| PREP-061 | tetracyanoplatinate(II) | C4N4Pt / -2 | 名称查询无匹配，保留 | d8 与线性配体 |
| PREP-062 | bis(trimethylphosphine)silver(I) | C6H18AgP2 / 1 | 名称查询无匹配，保留 | 线性局部环境 |
| PREP-063 | chlorotris(trimethylphosphine)copper(I) | C9H27ClCuP3 / 0 | 名称查询无匹配，保留 | d10 及 P/Cl 混合 |
| PREP-064 | hexamethyltungsten | C6H18W / 0 | 数据库字段不符，保留反例 | 实验畸变三棱柱候选 |
| PREP-065 | Co@Ge10(3-) | CoGe10 / -3 | 按原始结构来源继续核对 | 实验簇/环境核对；不预设理想 PPR-10 |
| PREP-066 | Pd@Bi10(4+) | Bi10Pd / 4 | 按原始结构来源继续核对 | 实验盐环境核对；不预设理想 PAPR-10 |
| PREP-067 | La(NO3)2(phen)(H2O)4(+) | C12H16LaN4O10 / 1 | 按原始结构来源继续核对 | 完整氢键环境与 TD-10 候选重新分类 |

源码数值/布局正在 NUM-FIX-005 固定版本复验，此目录只增加外部准备资料，不改变该轮输入、算法或判据。REF-001 保持用户暂停状态。
