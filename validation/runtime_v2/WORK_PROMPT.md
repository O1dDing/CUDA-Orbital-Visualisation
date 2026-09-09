# 给 Work 的交接 Prompt

继续维护 O1dDing/CUDA-Orbital-Visualisation 的 test/fchk-validation-native 分支。
先读取远端最新 HEAD 及 validation/runtime_v2/README.md、policy.py、resume.py、windows_job.py、test_runtime.py。
本次新增的是 REF-001 Runtime v2 执行层，不是 COV 渲染器修复，也不是正式科学验收。
请将以下规则写入本 Work 的持久项目记忆/续接记录，并在后续计算中遵守。

新的资源规则：根据 Windows 实际物理核心拓扑分配，SMT 线程不算额外物理核。
16 核及以上的机器合计最多给 Gaussian 14 核；少于 16 核按我的 98x3 小主机档位处理，
但仍以检测结果计算预算：通常 8 核机用最多 6 核、12 核机最多 10 核，保留系统余量，
不得因为 Windows 报告 16 线程而把 8 核机当成 16 核机。
按当前 REF-001 冻结基组的真实收缩壳层/NBasis、开闭壳层与空闲预算分配任务，最多两个并发，
总核数不能超过物理核预算。不能声称这些启发式已经是最优性能参数。
核数只在启动新子阶段/恢复 attempt 时变化，不热改运行中的 Gaussian；
保留 Gaussian 16W A.03 完整 startup affinity 的兼容处理，不重新引入 PGI 初始化崩溃。
每任务输入 24 GiB、进程树上限 32 GiB，并执行内存、磁盘、同机及旧版工作目录锁的安全检查。

恢复规则：Opt、Freq、Stable 分开持久化。有效检查点的未完成优化优先 Opt=Restart，保留实际 R/U 方法。
立即中断只终止本协调器拥有的 Job Object 进程树，等待所有后代退出再冷复制、formchk 和身份校验。
不能以 chk 存在、日志增长、SCF Delta-E 很小就认定恢复有效或几何已收敛。
解析频率不能仅靠 chk 恢复到任意 Hessian 中间进度；本版保留已完成优化，只重算解析频率，
不自动使用 Freq=Restart 或未经验证的 RWF Restart。稳定性/SP/Force 可从检查点 Guess=Read 重做该阶段，
但不宣称内部迭代零损失。坏检查点、重启历史拒绝、正常退出但缺关键证据必须 needs_review，保留原件。
没有明确科学通过时，candidate_collected 仍不等于 scientific_pass。

新版操作：4 是子阶段边界暂停、6 是撤销暂停并补满空闲 worker、8 是立即中断并校验检查点、
9 是收尾后退出协调器、3 是开始队列。暂停后协调器默认保留，不再出现一个 worker 永久退出无法补回。
旧版菜单没有被热更新，这些按键语义不能直接套到仍运行的 v1。

保留 validation/paused-20260906 全部冻结源码、输入、已有身份和证据，不覆盖原 work/jobs/REF-001。
v2 记录在 work/runtime-v2，已完成旧阶段通过冷导入复用，原始失败和中断 attempts 不删除。
模型方法、基组、电荷、自旋、SuperFineGrid、VeryTight、固定几何用途与两次 stability repair 上限不降低。
本次用户新授权只更新执行资源/恢复策略；不擅自放松科学协议。
真实 OLD-018 是 [Cr(NH3)6]3+，本次参考 NBasis 为 502，不是此前从 OLD-015 误带入的 422。

先做只读核对并运行测试，报告 GitHub HEAD 与测试结论。Python 模拟测试、全部冻结输入检查和
Windows Python Job Object 测试，都不得冒充真实 Gaussian 原生中断恢复验收。
接着在独立临时目录用小案例做真实 G16W 验收：优化中断/重启、核数改变、频率中断后仅频率重算、
R/U保持、双作业 pause/resume/interrupt、无外部进程误杀、CPU/内存预算和完整/恢复结果对照。
验收通过才部署到正式计算入口；先备份并检查本地未提交改动，不强制覆盖已有 Work 修改。
若 OLD-018 或旧协调器仍在运行，不擅自强杀、不删除 supervisor.lock、不自动重新计算全库；
先保存当前进程/输入/日志/检查点身份清单，明确迁移方案与实际验收情况。
只能在旧协调器和 Gaussian 全部退出并取得旧版同一路径 OS 锁后 import-legacy。
完成后记录新的启动路径、数据路径、资源策略、测试证据、待办与恢复边界，供下一次 Work 直接接续。
不要发 release、不要合并 main，也不要把尚未运行的原生验收写成已通过。
