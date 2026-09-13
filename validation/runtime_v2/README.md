# REF-001 自适应物理核运行时 3.0

本地 COV Helper 3.0 使用一套运行代码、配置、控制状态、协调器锁和持久资源账本。本实现随主程序开发分支维护，计算数据仍保持独立。

## 资源分配

| 主机物理核 | 总预算 | NBasis < 300 | 300–449 | ≥450 | 同时任务 |
|---|---:|---:|---:|---:|---:|
| 16 | 14 | 4 | 5 | 7 | 最多2 |
| 8 | 7 | 2 | 3 | 4 | 最多2 |

NBasis 来自冻结 REF 基组的收缩壳计数，并与已有同方法 Gaussian 日志交叉校验。所有任务保留精确交换、基组、电子态、VeryTight、SuperFineGrid 及原有 Opt/CalcFC/Freq 设置。档位是保守调度策略，不声称已测得性能最优。

整档原子预留在提交工作线程之前完成，不能按空余量降档。16核允许7+7；8核允许4+3、拒绝4+4。完成阶段排到等待队列末尾，无法容纳的队首不会被后续小任务绕过。RAM暂停仍占用原配额和一个任务名额；恢复同一进程无需释放再申请。当前阶段的 formchk 等串行工具复用同一预留；独立探针单独申请名额。

资源账本位于 `%LOCALAPPDATA%/COV/adaptive-physical-resources.json`，记录 requested/allocated/budget、并发数、策略版本、参考身份、PID及创建时间、所属Job、分配/归还/暂停事件。进程创建时由Windows直接纳入Job，避免协调器在创建和纳管之间崩溃留下孤儿。退出必须确认整个Job；崩溃后不凭超时清除仍活着的树。

物理核数取自Windows硬件拓扑，不把32个SMT逻辑处理器视为32核。每棵树还有Gaussian线程环境与Windows Job CPU时间硬上限。为兼容本机Gaussian16W A.03，仍保持完整启动亲和性；这不是固定物理核独占绑定，也不是连续占满7核的承诺。其他主机最多14核，2–8物理核保留1核、9核以上保留2核，并按预算比例向上取整缩放档位；多处理器组明确拒绝运行。

## 本地入口与目录

- Resume 是唯一部署目录；START.cmd 打开 COV-Helper.pyw 的桌面窗口。
- 每次打开只读取状态，手动计算控制默认关闭；刷新、部署检查及关闭窗口均不启动或恢复计算。
- 原独立部署的冻结计算材料迁入 Resume/recovery/frozen-reference-20260910，逐文件核验后旧目录移入回收站。
- Resume/runtime/active.json 指定当前不可变代码包；runtime/config.json 是共享配置。
- 新状态使用独立目录；实际目录由配置 runtime_directory 指定。旧 work/runtime-v2-fastpause、work/jobs/REF-001、原始归档、OLD-018-SALVAGE 保持原位。
- 新版部署和清理回执保存在 runtime/deployments；记录迁移材料、替换文件及回收站位置。
- 桌面操作、配置和目录说明见 [Helper 使用说明](HELPER.md)。下表仍列出源码命令行控制接口。

## 操作

| 命令 / 按键 | 行为 |
|---|---|
| check / 1 | 校验全部冻结输入和本机资源，不计算 |
| import-runtime / 2 | 在持有同机、同 work 锁时导入闲置的 2.1／2.2／3.0 结果和恢复来源 |
| run / 3 | 明确启动所选未完成案例，最多两槽，按物理核预算分配 |
| hold / 4 | 暂停所属计算线程，RAM 仍占用；不能关机或关闭协调器 |
| status / 5 | 查看状态、实际活动耗时、暂停诊断、最新收敛表 |
| resume / 6 | 原进程恢复，继续同一 attempt；超时后显式恢复获得下一段有界活动预算 |
| interrupt / 8 | 结束所属树，确认所有后代退出后冷保存及校验 CHK/RWF |
| shutdown / 9 | 恢复 RAM 暂停并等待当前子阶段收尾，随后退出协调器 |
| pause / 10 | 当前子阶段或已验收的优化检查点段收尾后暂停派发 |

暂停确认覆盖所属 Job 中的计算进程。Windows System32/conhost.exe 是系统控制台基础设施，按完整路径识别并保持运行；诊断列出该例外。未知控制台宿主使暂停失败并回滚。退出与资源限制仍覆盖整个所属 Job。

## 数据迁移与验收

import-runtime 先完整核对旧 binding、冻结参考身份、Gaussian 二进制和各文件 SHA，再发布导入元数据。完成阶段保留原生产者身份和原文件路径；开始下一阶段时独立复制检查点。中断来源额外绑定旧输入、manifest 和逐文件校验，任何更改会拒绝恢复。重复导入会校验并返回同一回执。迁移不启动 Gaussian，也不迁移旧 native-capabilities 回执。

原生验收只在独立的水 / NO 副本上运行。native-acceptance --capability ram_pause 只检验 RAM 暂停、恢复及最终结果与基线的一致性；不启动正式候选队列。完整 native-acceptance 还包含有独立门槛的 L103 分段和解析频率 RWF 恢复。没有对应运行时及 R/U 方法验收时，不能开启这些恢复能力。

## 开发、安装与回滚

保留 runtime_v2 与 paused-20260906 同级布局。独立源码运行时使用 --config 指定配置。

```text
python -B -X utf8 owned_test_runner.py --output runtime-test-results
python -B -X utf8 resume.py --config config.json check
python -B -X utf8 resume.py --config config.json import-runtime
python -B -X utf8 resume.py --config config.json native-acceptance --capability ram_pause
python -B -X utf8 install_runtime.py prepare --home <Resume> --config <配置>
python -B -X utf8 helper_deployment.py plan --home <Resume> --fast <待退役目录> --plan <计划文件>
python -B -X utf8 helper_deployment.py apply --home <Resume> --fast <待退役目录> --plan <计划文件> --receipt <部署回执>
```

Helper 部署要求运行时身份与已验收身份一致、RPBE1PBE 和 UPBE1PBE 的 RAM 验收均有效，且没有 Gaussian、协调器或已打开的 Helper。清理计划先固定文件身份，迁移并核验冻结材料，再退役旧目录。恢复旧文件前同样应关闭 Helper 并确认计算闲置；不要覆盖后续修改。所有计算记录和快照保留。

单任务输入24GB、Job树32GiB，Opt/Freq/Stable分开；candidate_collected仍不等于科学验收通过。273例完整分类及历史审计见[本次交付证据](../adaptive-scheduling-20260913/README.md)。原暂停机制背景见[暂停和恢复说明](FAST_PAUSE.md)及[日常操作](WORK_PROMPT.md)。
