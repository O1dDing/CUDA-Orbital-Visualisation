# 一次提交准备 COV 计算文件

[English](NBO_ONE_JOB.md)

[单次作业模板](../examples/nbo-one-job/cov_one_job.py)依次运行 Gaussian/NBO、`formchk` 和 GenNBO，将波函数和 NBO 文件整理到一个 COV 可打开的文件夹。GenNBO 分析已有波函数档案，不再进行一次 SCF。

模板从 COV 源码树运行，适用于 Windows、Python 3.11 或更新版本、Gaussian 16W A.03 和 NBO7 i8。Gaussian 与 NBO 需已安装，`gaunbo6.bat` 需已指向相应 NBO 安装目录。A.03 使用 `Pop=NBO6Read` 调用 NBO7；其他 Gaussian 修订版使用各自接口，本模板检查 A.03。

现有 v0.3.0 和 v0.4.0-pre.1 下载包不含这个模板。它使用源码树中保留的[进程管理文件](../tests/validation_process.py)。

## 提交计算

复制 [water.json](../examples/nbo-one-job/water.json)，填写分子、几何、方法、基组、电荷和多重度。示例为固定几何水分子，PBE1PBE / 6-31G(d)，电荷 0，多重度 1。需要自定义基组或 ECP 时，在 `basis_ecp_tail` 中填写完整的输入段。

在源码根目录运行，使用本机安装路径和一个尚不存在的输出目录：

```powershell
python examples/nbo-one-job/cov_one_job.py `
  --recipe examples/nbo-one-job/water.json `
  --output work/water-job `
  --gaussian-bin 'C:/Gaussian/Gaussian 16 W' `
  --nbo-bin 'C:/NBO/bin'
```

这条命令只预览输入、检查路径，不创建文件或开始计算。在同一命令后加 `--run` 即提交作业。模板使用给定几何进行单点计算；优化和频率计算需另行准备输入。

最多使用三个 Gaussian 线程，整棵进程树限制为三个核心、20 GiB。Gaussian、`formchk`、GenNBO 的时限分别为 900、300、600 秒。`--memory-limit-gib` 和 `--gaussian-timeout` 可调整相应限制；`--help` 列出选项。

## 打开结果

在 NBO 预览版中打开 `work/water-job/cov-package`，或将其中的 `drop.covnbopkg` 拖入 COV。该目录包含 `canonical.fchk`、一份 `analysis.nbo` 报告、`.47` 档案和请求的轨道矩阵。

原始 canonical checkpoint 和 Gaussian 档案保留。输入目录中的报告、档案和矩阵来自同一次 GenNBO 分析；运行记录留在目录外。失败时保留输出，再次尝试使用新目录。

## 各视图的数据

FCHK 和 NBO 数据应对应同一步计算。档案提供重叠、密度、canonical MO 和 Fock 数据，报告提供轨道身份和打印的数值。下列编号是模板设置，COV 按矩阵标题读取。

| 想查看的内容 | 模板生成的矩阵数据 |
| --- | --- |
| MO 能量、占据、Gaussian AO–MO 系数联系和三维 MO | `canonical.fchk` |
| NAO 形状、NPA 电荷和 Wiberg 着色 | AONAO（W52）、报告及对应密度 |
| MO 的 NAO 组成 | AONAO（W52）、NAOMO（W51）；NAONBO（W53）提供额外变换 |
| NBO 形状及 MO–NBO 联系 | AONBO（W37）、NBOMO（W49） |
| NHO 形状和组成 | AONHO（W54）、NAONHO（W57）、AONAO（W52）；NHONBO（W59）用于连接 NBO |
| NLMO 形状、主要 NBO 和尾部 | AONLMO（W55）、NBONLMO（W60）、AONBO（W37）；NLMOMO（W61）用于连接 MO |
| PNAO 形状 | AOPNAO（W56）、AONAO（W52） |
| E(2) 供体–受体轨道 | AONBO（W37）、报告和档案中的 Fock 数据 |
| NAO/片段能量和 SALC 数据 | AONAO（W52）、对应密度和 Fock 数据；可用对称性取决于几何与所选片段 |

模板另写出 SAO、NAONLMO 和 AOMO。缺少数据时，相应视图不可用；开壳层自旋显示需要实际的分自旋数据。

矩阵及档案选项见 [NBO 7 手册](https://nbo.chem.wisc.edu/nboman.pdf)。COV 的视图和导出见[使用 NBO 结果](AOMO_NBO.zh-CN.md)。
