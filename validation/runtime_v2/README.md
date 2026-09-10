# REF-001 Runtime 2.1：快速暂停与检查点恢复

当前操作规范以 [FAST_PAUSE.md](FAST_PAUSE.md) 为准。旧版 README 中的 4 键和
`work/runtime-v2` 路径属于历史 v2；本版运行证据独立写到 `work/runtime-v2-fastpause`。
不热升级旧进程，不覆盖 `paused-20260906` 或旧计算结果，不发布 COV release。

## 快速开始

在完整源码目录运行，保留 `runtime_v2` 与 `paused-20260906/runner-source` 的同级关系。
复制 `config.example.json` 为 `config.json`，核对现有 data/work/Gaussian 路径。
双击 `START.cmd` 或用原机 Python：

```powershell
F:\Dev\Python312\python.exe -X utf8 resume.py check
F:\Dev\Python312\python.exe -X utf8 resume.py native-acceptance
```

第二条会启动独立临时小案例的真实 Gaussian 验收，不使用 OLD-018 抢救检查点做试验。
它记录 RAM pause、Opt 分段和解析频率 RWF restart 的实际能力；不是模拟测试。
本地 Gaussian 安装/方法没有对应验收回执时，不自动启用 KJob 分段或 RWF 恢复。
CI 使用 Python 子进程及模拟数据，**不能替代这一步原生验收**。

旧 v1 退出并通过同路径 OS 锁检查后可 `import-legacy` 冷导入有效完成阶段。
已有旧 v2 运行结果或 OLD-018 的历史 timeout，需要 Work 做明确的只读核对、
reviewed migration/salvage 导入；不得只解除 needs_review 就从初始几何重算。

## 按键

| 按键 | 行为 |
|---|---|
| 3 | 开始队列，自动分配核心 |
| 4 | **快速内存暂停**，不等待整段 Opt/Freq；内存仍占用，不可关机或关闭协调器 |
| 6 | 原进程继续；恢复派发并补上空闲计算槽位 |
| 10 | 已验收的优化检查点段／当前子阶段收尾后暂停 |
| 8 | 输入 SAVE 确认后，停止所属树并冷保存检查点；未写入的工作可能损失 |
| 9 | 当前段收尾后退出协调器；RAM 暂停时会先恢复以完成边界 |
| 5 | 看活动树、暂停状态、真实耗时、最近三组几何收敛表 |

RAM 暂停不刷写应用缓存，不是磁盘保存。8 的复制、刷盘、哈希和 formchk 完成前不能
宣告保存完毕；保存的 bytes 可读也不等于 Gaussian 优化历史／RWF 已证明能恢复。

## 资源与科学协议

按 Windows 实际物理核心拓扑分配，不把 SMT 当额外物理核。16 核以上总预算最多14；
12核为10，8核为6，6核为4，4核为2，2/1核为1。少于16核标为用户的98x3小主机档位，
这只是档位名称，不伪称识别了实际 CPU 型号。多处理器组需要单独审查，当前明确拒绝。

单任务 1–14 核，以本轮真实基组 NBasis、开闭壳层、空闲预算分配；最多两个并发。
每任务输入24GB、Windows Job树32GiB，至少48GiB可见内存；派发还受可用内存和磁盘余量约束。
默认72小时计入实际活动墙钟时间，RAM暂停不计入；期限到达时默认hold，不直接杀Gaussian。
Gaussian 16W A.03完整启动affinity兼容策略保留，预算不是物理核独占绑定。

Opt/Freq/Stable分开，必要时最多两次稳定性修复。基组、模型、电子态、VeryTight、SuperFineGrid、
固定几何用途保持不变。原子是SP，固定几何／扫描是Force。candidate_collected不等于科学验收。
Opt重启优先有效CHK；解析Freq细粒度恢复用同生产者整套命名RWF，不能用Freq=Restart替代。
详见 FAST_PAUSE.md 的边界、实际验收门槛和技术来源。

## 回归与交接

```powershell
F:\Dev\Python312\python.exe -X utf8 test_runtime.py
F:\Dev\Python312\python.exe -X utf8 test_fast_pause.py
```

给 Work 的持久接续指令见 [WORK_PROMPT.md](WORK_PROMPT.md)。所有 attempt 和失败证据保留，
控制命令有单独日志，不删除锁绕过互斥，不按进程名批量杀 Gaussian，不静默降精度或无限重试。
