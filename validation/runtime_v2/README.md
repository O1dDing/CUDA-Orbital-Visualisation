# REF-001 统一运行时 2.2

两个本地入口共用一套运行代码、配置、控制状态和协调器锁。正式程序 PRE11 与本验证运行时分别发布。

## 本地入口与目录

- Resume/START.cmd 保留原入口；cov_resume_menu.py 兼容旧 --execute、--limit、--workers、--case 参数。
- FastPause-2.1/validation/runtime_v2/START.cmd 和 resume.py 转到同一版本。
- Resume/runtime/active.json 指定当前不可变代码包；runtime/config.json 是共享配置。
- 新状态使用独立目录；实际目录由配置 runtime_directory 指定。旧 work/runtime-v2-fastpause、work/jobs/REF-001、原始归档、OLD-018-SALVAGE 保持原位。
- 被替换的入口和代码保存在 runtime/backups，每份替换回执记录原件、备份及新文件的 SHA-256。

## 操作

| 命令 / 按键 | 行为 |
|---|---|
| check / 1 | 校验全部冻结输入和本机资源，不计算 |
| import-runtime / 2 | 在持有同机、同 work 锁时导入闲置的 2.1 结果和恢复来源 |
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
python -B -X utf8 install_runtime.py activate --home <Resume> --fast <FastPause-2.1> --config <配置>
python -B -X utf8 install_runtime.py rollback --home <Resume> --config <配置> --receipt <替换回执>
```

激活要求新运行时已导入数据、RPBE1PBE 和 UPBE1PBE 的 RAM 验收均通过，且当前没有 Gaussian/协调器。回滚同样要求闲置；若入口在激活后被修改则拒绝覆盖。所有计算记录和快照保留。

资源与科学协议保持原值：物理核总预算最多 14，单任务输入 24GB、Job 树 32GiB，Opt/Freq/Stable 分开，既有基组、模型、电子态、VeryTight 和 SuperFineGrid 不变。candidate_collected 仍不等于科学验收通过。细节见 [暂停和恢复说明](FAST_PAUSE.md)、[日常操作](WORK_PROMPT.md) 和 [修复证据](../runtime-integration-20260913/README.md)。
