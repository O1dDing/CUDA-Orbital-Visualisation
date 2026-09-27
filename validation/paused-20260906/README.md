# COV 原库与 REF-001 暂停交接（2026-09-06）

本目录和配套 Release 保存实际输入、实际进度及恢复工具。这是续算资料包，不是通过验收的软件版本。

2026-09-06 用户调整执行顺序：暂停 REF-001，日后分批续算；当前先实施所有不依赖质量参考完成的已确认通解修复及适用验证。旧方案、契约中“先完成全部 REF 再改生产算法”的执行门槛由此取代，历史文件本身保留。几何、电子态、环境和竞争方法等依赖新参考的结论仍待验证。

| 项目 | 暂停时实际状态 |
|---|---|
| 固定原库 | OLD-001 至 OLD-273，编号不变 |
| 原方法独立复现 | 273 份全部采集；不等同全部科学通过 |
| REF-001 | 13 个候选采集完成，2 个中断保留，258 个未开始 |
| 中断位置 | OLD-013 的 relax-after-stability-01；OLD-015 的 initial |
| 正式完整通过 | 0；不能将暂停、未实现、缺证记为通过 |
| 正式外部分子 | 0；≥50 个外部分子的目标尚未完成 |
| 当前计算 | REF-001 已停止；已确认协调器及作业退出，资源锁可重新取得 |

`case-status.json` 是暂停时逐案例的实际状态。原始 `stage.json` / `result.json` 保存协调器最后一次写入，部分仍显示历史 `running`；它们不证明现在有进程运行。恢复工具会按照原有规则将这些尝试记为中断，使用新的尝试目录。

源码分支：[test/fchk-validation-native](https://github.com/O1dDing/CUDA-Orbital-Visualisation/tree/test/fchk-validation-native)。
数据发布：[cov-ref-paused-20260906](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906)。

## 包含的实际材料

- `reference-inputs/OLD-*/initial.gjf`、`basis.gbs`、`reference.json`：全部 273 个已准备参考输入，含冻结的基组/ECP。
- `runner-source/`：REF-001 实际使用的四个源码文件及源码身份。生产主分支里的新调度修复不会偷偷改变这份历史运行身份。
- `campaign.json`、`reference-candidates.json`、`case-status.json`、问题台账和通解契约：案例身份、物理用途、裁决和待办。
- 23 个 ZIP 数据分卷及 `data-manifest.json`：原始 273 FCHK、既有原库原生取证、273 独立复现、独立数值参考、REF 已完成和中断的所有检查点、临时数据与日志。逐文件及逐分卷均有 SHA-256。
- `resume_reference.py`：校验/解包、迁移工作路径、少量续算、阶段边界暂停。默认 `run` 只选择一个未完成案例。

Gaussian 软件本体及第三方运行库、编译中间件、缓存和失效锁不在资料包内。需使用自己已安装且有许可的 Gaussian 16W；同一 REF-001 要求 `g16.exe` 和 `formchk.exe` 的内容身份与原运行一致。更换 Gaussian、输入或科学策略时应建立新轮次，保留旧结果，不绕过身份检查。由计算生成的结果、检查点和残留临时文件均已纳入数据清单。

## 校验与恢复

使用 Python 3.12。校验、解包和准备不需要运行 Gaussian。将 Release 中全部 ZIP 分卷下载到同一 `assets` 目录，并取得本目录的 `data-manifest.json`。

```powershell
python .\resume_reference.py verify --manifest .\data-manifest.json --assets D:\COV-archive\assets --data-root D:\COV-archive\data --extract
python .\resume_reference.py prepare --data-root D:\COV-archive\data --work-root D:\COV-incremental
python .\resume_reference.py check --data-root D:\COV-archive\data --work-root D:\COV-incremental --gaussian "F:\CalChem\Gaussian\Gaussian 16 W"
```

第一步核验每个分卷及其每个文件，安全解包到相对布局。第二步复制 REF 工作材料，重定位工作状态中的路径。第三步只校验源码、Gaussian、273 个实际输入和所有已采集阶段的身份，不启动计算。归档证据不改写；路径不是只有历史机器上的绝对引用。不要将数据目录和续算工作目录设为相同目录或彼此的子目录。

当前暂停保持有效；以下命令仅在以后明确决定续算时使用：

```powershell
python .\resume_reference.py run --data-root D:\COV-archive\data --work-root D:\COV-incremental --gaussian "F:\CalChem\Gaussian\Gaussian 16 W" --limit 1 --workers 1
```

也可以用 `--case OLD-013` 指定未完成案例，或明确增加 `--limit`。一次调用的选择清单写入 `incremental-runs/<运行编号>/progress.json`。已经采集的案例不重复计算；已记录的失败需要分析后建立经过审阅的新尝试/轮次，不能默默重试直到成功。

正常暂停采用另一终端发出：

```powershell
python .\resume_reference.py pause --work-root D:\COV-incremental
```

正在执行的阶段会收尾，包括适用的 FCHK 转换；新 Gaussian 阶段和新案例不再派发。查看运行状态并确认进程已退出后，未来要再续算时明确移除工作目录内的 `PAUSE.json`。强制终止或断电留下的尝试不会当作完成结果；默认从已核验的上一完整阶段检查点重新运行中断阶段，初始阶段则重新使用独立初始猜测。部分写入的检查点和临时文件作为中断证据保留，不自动晋升为参考。

## 资源与判据

每个 Gaussian 输入为 4 核、24 GiB；整个进程树内存硬上限 32 GiB，CPU 硬配额对应 4 个逻辑处理器。最多 2 个并发作业，为协调、分析或编译保留余量。全部本任务计算始终合计不超过 12 核、128 GiB；不要在其他满额编译/验证同时开启参考作业。Gaussian 16W 的已验证运行方式保留完整处理器亲和性，以输入线程数和 Windows Job 配额限额。

保留冻结的 PBE0-D3(BJ)、def2-TZVPP/TZVPPD、严格 SCF、细网格、优化/频率/稳定性及预设收敛重试规则。一次正常结束或局部优化不代表真实全局基态。REF-001 全库物理裁决尚未完成，竞争电子态、环境及适用交叉方法仍待补足。

当前源码修复按共同根因实施，原子顺序、分子名和固定 MO 号不作特判。每个生产版本冻结后对完整 273 取证，全部进入终止状态后统一分析；外部集合加入后每轮覆盖全正式集合。实际纹理全 MO 同点门槛仍为 NRMS≤1e-4、绝对余弦≥1−1e-7、相对最大误差≤1e-3。真实界面交互、导出一致性、独立积分、密度、数值不变量及普通版/验证版一致性均保留原验收要求。

## 当前修复顺序

1. AO 相位/顺序和 CPU 精度，明确 GPU 上传精度；基组独立解析重叠积分与矩形 MO。
2. 分解累加初值、可用状态和作用范围；保留合法可读性规则，验证守恒及敏感性。
3. 整体/局部/近似对称性、简并子空间、拓扑与次级联系。
4. 成键/配位场/π 范围、实际控件、文字与屏幕及导出快照。
5. 每版完整 273 取证和统一交叉分析；物理参考缺项保持待证，不阻塞以上独立修复。

只有全部正式案例的适用检查及所需证据齐备，并交付同源码普通版和验证版后，才能宣称完整验收完成。本暂停资料包不做这种声明。
