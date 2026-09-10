# REF-001 快速暂停版 2.1 — 2026-09-10

本文是当前快速暂停入口的操作说明；原 README 的 4 键和 work/runtime-v2 路径只适用于旧版。
新代码仍位于 validation/runtime_v2，运行数据另写 **work/runtime-v2-fastpause**，不会冒充或覆盖旧 v1/v2 结果。
只改验证分支；不热接管旧 Gaussian、不发布 COV release、不改冻结科学输入。

## 核心区别：暂停不等于保存

| 按键 / 命令 | 行为 | 是否能关闭计算窗口/关机 |
|---|---|---|
| **4 / hold** | 快速暂停所属 Job 内所有线程，保留 RAM 中的计算现场；不等待整段 Opt/Freq | **不能**。不是磁盘检查点；内存仍占用 |
| **6 / resume** | 原进程原地恢复，继续同一 attempt；填回空闲计算槽位 | 仍在计算 |
| **10 / pause** | 在已启用的优化检查点段或当前子阶段结束后暂停派发 | 等状态确认没有活动树、文件校验完成 |
| **8 / interrupt** | 明确确认 SAVE 后结束本协调器拥有的树，再冷复制/校验 CHK；Freq 同时保存 RWF/INT/D2E | 等冷保存和 probe 回执完成；恢复仍有明确边界 |
| **9 / shutdown** | 当前子阶段/检查点段收尾后关闭协调器；若在 RAM 暂停，先恢复计算去完成该边界 | 等协调器退出 |
| **5 / status** | 看真实 execution_state、暂停时长、活动耗时、线程句柄、几何收敛表 | 以真实状态为准 |

**RAM 暂停是本次解决“点暂停后等数小时”的默认快捷路径。**系统调用轮询周期为 0.1 秒，
正常树的确认需要数次短扫描，但并不承诺任何主机固定 0.1 秒完成；具体延迟要看 Windows 原生测试。
这不是程序内的事务式保存，内核 I/O 仍可能收尾。**绝不在 RAM 暂停时复制 CHK/RWF 或运行 formchk**。
保留内存、已分配核数预算和 Job 句柄，不能用暂停期间的空闲 CPU 再超额启动一批任务。
暂停/恢复只作用于持有的 Windows Job 及其后代；不依据进程名、陈旧 PID 或全局 taskkill。
不使用未经文档化的 NtSuspendProcess/Job freeze API。

微软明确提醒 SuspendThread 不是线程同步工具。这里是独立外部监督器的 debugger-style 控制，
不获取目标程序拥有的 mutex/临界区，不要求被暂停线程为监督器完成锁操作。
逐线程持有真实句柄、核对所属 Job、记录唯一一次 suspend increment，恢复时只撤销自己那一次。
暂停扫描不能稳定时撤销部分暂停并报错，不偷偷把暂停失败改成杀进程；所有路径仍需原机 Gaussian 验收。

## 超时不再把长任务直接砍掉

受管 Gaussian 的 72 小时预算按 **active elapsed time** 计算，RAM 暂停时长不计入。
预算到期默认请求 RAM hold，而不是 TerminateJobObject；随后由操作者选择继续、边界停止或冷保存。
界面记录 wall_seconds、active_wall_seconds、ram_paused_seconds。formchk 等工具仍有独立有界超时。
RAM hold 仍须保持计算窗口、系统电源和足够内存。断电不能靠这个模式恢复。

## 快速落盘与 RWF 恢复

新 attempt 使用明确的 `%RWF=job.rwf`、`%Int=job.int`、`%D2E=job.d2e` 和 `%Save`。
其路径属于每个独立 attempt，不复用上一 attempt 的可写文件。原件不自动清除；磁盘空间须留余量。

8 键的顺序是：停止所属进程树 → 等全部后代退出 → 冷复制 → fsync/哈希 → 原子提交 manifest → CHK probe。
CTRL+C 仍按显式中断处理，不等同于 RAM hold。命令保留 interrupt_epoch：紧接着按 6 不能让旧 attempt 漏掉已请求的中断。

优化通常只复制 CHK、输入和日志以降低保存时间，原 attempt 的 RWF 保留；频率则保存命名的 RWF/INT/D2E/CHK 整套。
快照标记 `restart_validated=false`：**文件保存正确不等于 Gaussian 已证明能重启**。
磁盘不够、哈希改变、CHK 不可读时原件保留，明确报错；不输出虚假的“安全保存完成”。
大 RWF 的冷复制、刷盘、校验可能花几分钟；保证快速释放计算的机制是 RAM hold，不能承诺任意大小文件瞬间落盘。

解析频率从有原生验收回执的快照复制到新 attempt 后，生成 `#p Restart`。
不生成 `Freq=Restart`（那是数值频率恢复路径），不把文件搬到别处后未经检查直接续算。
每次核对模型/基组/电荷/电子态、runtime identity、生产者二进制指纹和逐文件哈希。
Gaussian 安装中的 link EXE/DLL 也进入新身份，避免只换 l1002 而 launcher 不变被误认为同一版本。
RWF 被拒绝或没有验收证据时进入 needs_review；只有显式选择 `freq_recovery="replay"` 才允许仅重做频率。

## 优化步级磁盘检查点：有原生门槛，不能仅凭日志行热杀

对通过本机实测的安装/方法，`opt_step_checkpoints="auto"` 会在新的优化 attempt 中预置 `%KJob L103 2`。
这是让 Gaussian 自己在第二次 l103 返回后终止当前执行，再校验/冷保存/Opt=Restart；不是检测到
“Maximum Force”后抢时间杀进程，也不修改正在运行的 .gjf。

只有实际日志满足严格的两次 l103 返回、优化证据且无普通 SCF/IO 错误，才记录 `checkpointed`。
它不是 collected/candidate_collected。Opt 完成另有显式完成证据；模型化学、VeryTight、网格不放松。
需验证该安装的 KJob 终止格式、实际优化历史前进和最终结果一致。遇到未知终止格式，明确拒绝。
这会带来更多启动/校验和磁盘开销；等待时间仍可能包含一个昂贵几何步，不承诺秒级落盘。
达到 1024 段仍不完成会进入 needs_review，防止无进展自动重启循环。

## 使用和原生验收

首次在新目录部署完整源码与冻结依赖，保留现有 data/work；不要覆盖正在运行的任何入口。

```powershell
F:\Dev\Python312\python.exe -X utf8 resume.py check
F:\Dev\Python312\python.exe -X utf8 resume.py native-acceptance
```

`native-acceptance` 必须取得同机锁、旧目录锁，确认没有其他 Gaussian；只在
`work/runtime-v2-fastpause/native-acceptance/<唯一目录>` 跑小案例水和 NO 的独立试验副本。
初始几何在副本中拉伸12%，保留本轮基组和数值精度，显式标为验收试验，**不算正式候选结果**。
测试完整 baseline、RAM 暂停继续、KJob 分段 Opt=Restart/换核数、解析频率中断/RWF Restart，
比较能量、原子间距离和频率。某项太快而未触发中断不算通过。
输出真实日志/哈希与 native-capabilities.json；只为实际通过的 R/U 方法开放对应功能。
任一失败不会伪造 native pass；KJob 失败不阻挡独立调查 RWF。
自动测试与 Windows Python 父子进程测试都**不是**这一步真实 Gaussian 验收。

`config.example.json` 使用 `opt_step_checkpoints="auto"`、`freq_recovery="auto"`：
没有 Opt 原生回执时不注入 KJob，仍能快速内存暂停；有回执才自动分段。
有 RWF 却无对应恢复回执时不静默丢弃计算转为全频率重算。强制 `true` 或 `rwf` 也不能绕过门槛。

旧 v1 通过 `import-legacy` 冷导入已完成阶段；新版独立目录不覆盖旧 runtime-v2。
旧 v2 已有运行成果或 OLD-018 这种旧 timeout/failed 的恢复，仍需 reviewed migration，
不能只删 review.json 再运行，避免把未导入的断点变成 fresh input。原 OLD-018-SALVAGE 不动。

新菜单 4 与旧菜单 4 不同，显示明确的“内存暂停，不可关机”。日常临时腾出算力用 **4 → 6**；
确认要退出时用 **10/9** 或了解损失边界后用 **8**，并查看保存回执。

## 来源（2026-09-10 检索）

Gaussian 主站访问 502，以下 Gaussian 手册镜像用于核对相同关键字定义：
- Gaussian Link 0 / %KJob, %Save, %RWF： https://www.conflex.co.jp/gaussian_support/link0.php
- Gaussian Restart（解析频率 RWF、不是 CHK-only）： https://www.conflex.co.jp/gaussian_support/restart.php
- Gaussian Freq： https://www.conflex.co.jp/gaussian_support/freq.php
- Gaussian SCF： https://conflex.co.jp/gaussian_support/scf.php
- Windows SuspendThread（同步/死锁限制）： https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-suspendthread
- Windows Job process ID list： https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-jobobject_basic_process_id_list
- Windows ResumeThread： https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-resumethread

恢复粒度不是统一的“每步绝对安全保证”。CHK 能被 formchk 读取，只证明可解析；
优化历史、同一模型化学的重启有效性和科学收敛仍须各自验证。
