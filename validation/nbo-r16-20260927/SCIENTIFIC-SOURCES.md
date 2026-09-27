# R16 科学语义依据与实际验证的对应关系

本轮使用已保存的量化输出，不新增计算。文献确定解释边界；具体数值由本地原始 FILE.47/.51/.52、FCHK 和报告独立复算，不能用文献中的另一套计算值替代。

| 第一方资料 | 对应实现与验收约束 |
|---|---|
| [NBO 官方 NPA 说明](https://nbo.chem.wisc.edu/ex_npa.htm) | NPA 总电荷由原子 NAO 占据求和得到。O₂ 的总电荷与 α−β 自旋分别检验，不允许自旋记录覆盖总电荷。 |
| [NBO 官方 MO/NBO 教程](https://nbo.chem.wisc.edu/tut_cmo.htm) | MO 可混合多种局域轨道；用 NBOMO/NLMOMO 等变换建立关系，不能由局域 NBO 的角色直接断言整个 canonical MO 的给受体角色。 |
| [NBO 官方 NLMO 说明](https://nbo.chem.wisc.edu/ex_nlmo.htm) | 母 NBO 与离域尾部是不同的可追踪组成。尾部不是自动成立的多中心键，需保留完整变换和相应证据状态。 |
| [NBO 官方理论说明](https://nbo6.chem.wisc.edu/webnbo_css.htm) | σ/π 相对于局部原子对轴定义；同一原子对可有多个轨道通道，不能把整条键压成一个文字标签。 |
| [NBO 官方一电子性质分析](https://nbo.chem.wisc.edu/ex_prop.htm) | Fock/KS 算符与总能量不是同一对象。本实现的交叉矩阵项不冒称 π 键能或严格总能量分解。 |
| [NBO7 官方手册](https://nbo.chem.wisc.edu/nboman.pdf) | 区分原始报告、矩阵变换、NPA、Wiberg 和 E2 的方法及来源；阈值未打印不等于物理零。 |
| [TAMU Hughbanks 群论及投影讲义](https://www.chem.tamu.edu/rgroup/hughbanks/courses/673/lecturenotes/notes-6_689.pdf) | 保留完整简并子空间；独立审查采用投影子空间、耦合奇异值及整个组的交叉算符特征值区间，避免任意伙伴旋转改变结论。 |

实际数值参考：`science/audit-results.json`、`cr-reference/reference.json` 与 `reference/*`。TiF₆ 检查集体金属 d–全部配体 pπ 空间，乙炔检查任意空间朝向的两端 π 空间，Cr(CO)₆ 由独立原始小矩阵与实际 E2 端点核对相反回馈方向；同时用 O₂、CO₂、XeF₂、H₂O 检查组成与能力路由。产品是否通过应读最终验收记录，不能由这份资料清单判断。
