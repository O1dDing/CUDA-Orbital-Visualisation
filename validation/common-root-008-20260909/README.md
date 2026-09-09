# COMMON-ROOT-008 稳定候选接续记录

生产程序和本轮采集源码固定在 `8e206e188ef84a5e63129ee20e5b95c2f5ad498b`。
四组共同根因修正已经完成适用专项控制和固定 40 例开发预检；流程门允许新全 273。
这份记录保存的是启动检查点，后续实际进度须读取下述原始轮次，不能把启动记录当终态。

## 已完成的证据

- [固定40统一裁决](common008c-precheck-adjudication.json)：整体通过，39123次界面动作零失败，7791个采样MO、178份完整前沿纹理、320组导出；数值、视图和声明模型范围各40例通过。采集1722.383秒；萘OLD-015采集及封装40.005秒。独立参考200次命中，当前COV输出均重新比较。
- [完整专项控制](common008c-complete-local-controls.json)：普通/验证构建各38项CTest，8个数学探针、原生视区、20张实际可用性/自旋/语言截图，以及两份g壳层完整实测。
- [g壳层裁决](common008c-g-control-adjudication.json)：2038次动作，500个采样MO、12份完整前沿纹理、16组导出通过；复用既有输入，未启动Gaussian。
- [固定40代表性实图](fixed40-visual-review.json)及[同源码20帧专项实图](availability-visual-review.json)：实际目检完成。未知、不适用和测得零值分开；alpha/beta字形、局部标签与整体来源分开记录。
- [流程门回执](common008c-workflow-gate.json)及[冻结前计划](candidate-plan-frozen.json)：当前完整控制、40项实测与实图身份核验通过。前两轮失败预检仍完整保留。

## 运行中的新全量轮次

[启动回执](common008c-original-start.json)对应 `COMMON-ROOT-008C-original`，全部273个原库文件。
轮次身份：`8c525e0e7421997dd5be928b6fea6c302c886fb98102cbe629e88865c78ed96c`。

原始轮次目录：`F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906\COMMON-ROOT-008C-original`。
该目录的 `progress.json`、`collection-complete.json` 和后续各复核目录是实际终态依据。
活动协调记录在 `F:\Codex\2026-09-09\new-chat\outputs\common-root-008-integration-20260909\common008c-original-review-chain`。

初始采集会话41874、复核会话18534、复核协调PID2976；恢复时先核实实际进程，不能仅凭历史PID操作。
全部采集终态后才依次数值、视图、声明模型范围复核；随后统一裁决。程序、判据与输入保持冻结。
仅COV总预算4物理核/64GiB，两个工作进程各2核/24GiB；没有启动或更改Gaussian队列。
清单中的历史独立计算预留字段不能替代用户后来确认的Runtime v2资源政策，也不授权启动Gaussian。

## 复核器修正与验收范围

g壳层的旧视图汇总仅为OLD文件名登记哈希，两个AUX报告因此未进入身份清单。
[修正身份记录](auxiliary-review-correction.json)保留原脚本和失败位置；新view-review-v2只复核原始证据，补齐两个实际案例的报告身份。
本提交把同一行修正同步到源复核器。逐案例检查函数与阈值没有变化；本次全量轮次的冻结文件不受影响。

当前通过仅覆盖已声明的数值、界面和数学模型范围。局部投影吻合不独立证明物理电子态、几何或化学机制。
完整科学通过仍为0。原273加至少50个不同外部分子的最终验收、物理参考裁决与最终EXE交付仍未完成。
旧NUM-FIX-007C不重跑；不覆盖旧失败或原始材料；不发布COV release，不更改正式标签或main/feat分支。

文件身份见[artifact-index.json](artifact-index.json)。这里只保存小型回执，全部原始证据保留原位置。
