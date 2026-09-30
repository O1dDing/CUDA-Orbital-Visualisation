# Chemical Orbital Visualiser (COV)

[English](README.md) · **简体中文** · [日本語](README.ja.md) · [Français](README.fr.md)

查看轨道能量、电子占据情况，以及原子轨道与分子轨道之间的联系。

COV 可以用交互式能级图查看轨道能量和电子占据，并导出图片与数据。软件支持 Gaussian FCHK/FCH 和 Molden 文件。

v0.4 预览版加入了 NBO 分析。配合计算生成的 NBO 文件，可以查看原子轨道和局域轨道怎样组成分子轨道，沿着图中的连线查看各项贡献，以及电荷和成键信息。

选中轨道后，还可以查看它的三维形状。轨道显示由 NVIDIA GPU 加速，界面支持中文、英文、日文和法文。

## 功能

- **查看能级图。** 查看轨道能量和电子占据情况、切换能量单位，并将能级图导出为 PNG 或 SVG。
- **追踪轨道联系 — v0.4 预览版。** 查看原子轨道和局域轨道怎样组成分子轨道，点击连线查看各项贡献。
- **查看电荷与成键 — v0.4 预览版。** 如果 NBO 文件中包含相关数据，可以查看 NPA 电荷、自旋布居、Wiberg 键级指数和给体–受体相互作用。
- **导出图像和数据。** 将能级图保存为 PNG 或 SVG，并将相应数据导出为 CSV 或 JSON。
- **浏览轨道。** 跳转到 HOMO 或 LUMO、搜索轨道列表，并显示已占据、未占据、内层或价层轨道。
- **查看三维轨道。** 选择轨道、调整其等值面，并旋转或缩放分子视图。

## 下载

| 版本 | 包含的功能 | Windows 下载 |
|---|---|---|
| [稳定版 v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | 轨道能量、占据情况、能级图和三维视图 | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [预览版 v0.4.0-pre.1](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.1) | 增加 NBO 分析、轨道组成和轨道联系 | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.1/CUDA-Orbital-Visualisation-v0.4.0-pre.1-Windows-sm120.zip) |

现有下载包仍显示旧产品名称。

现有 Windows 包适用于 RTX 50 系列；其他架构的安装包尚未提供。

## 开始使用

1. 下载并解压 Windows ZIP。v0.3.0：运行 `cov.exe`。v0.4.0-pre.1：运行 `program/cov.exe`，或双击 `Open-*.cmd` 打开样本。
2. 打开 Gaussian FCHK/FCH 或兼容的 Molden 文件，也可以将文件拖入窗口。
3. 在轨道列表和能级图中查看能量与占据情况。选择一个能级以查看对应轨道，或导出能级图。
4. 使用 v0.4 预览版时，可以打开计算文件夹，一并加载其中的波函数和 NBO 文件。如果文件夹中有多次计算，请选择要打开的一次。

已安装 `formchk` 时可直接打开 CHK；FCHK/FCH 可直接打开。NBO 文件与 FCHK 应来自同一步计算；采用 GenNBO 重新分析时，报告、档案和轨道矩阵应来自该次分析。

## 输入文件与运行要求

- **波函数：** Gaussian `.fchk` / `.fch`，或兼容的 `.molden` / `.mol` / `.input` 文件。已安装 `formchk` 时也可打开 Gaussian `.chk`。
- **NBO 文件 — v0.4 预览版：** 报告、归档文件和轨道矩阵分别提供分析的不同部分。各视图所需的文件详见[使用 NBO 结果](docs/AOMO_NBO.zh-CN.md)。
- **分子大小：** 每个输入最多包含 100 个原子。
- **图形环境：** NVIDIA GPU、兼容的驱动程序，以及 OpenGL 2.1 或更新版本。请使用适合显卡架构的程序版本。

## 文档

- [使用 COV](docs/UI.zh-CN.md)
- [使用 NBO 结果](docs/AOMO_NBO.zh-CN.md) — v0.4 预览版
- [一次提交准备计算文件](docs/NBO_ONE_JOB.zh-CN.md) — 源码树模板
- [源码编译](docs/BUILD.zh-CN.md)
- [稳定版发布说明](docs/releases/v0.3.0.md) · [预览版发布说明](docs/releases/v0.4.0-pre.1.md)

[反馈问题或建议功能](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose)。

## 许可证

[Apache License 2.0](LICENSE)。随附程序库的许可证列于[第三方声明](THIRD_PARTY_NOTICES.md)。
