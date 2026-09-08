# REF-001 Runtime v2：中断恢复与物理核预算

## 范围与状态

这是 `test/fchk-validation-native` 的新增执行层，不是 COV 渲染器改动，也不是新 release。
`validation/paused-20260906/` 的冻结源码、输入、已有科学身份及原始证据保持原字节不变。
v2 的来源身份单独记录在 `work/runtime-v2/binding.json`，运行结果不冒充旧 runner 产生的结果。
本执行策略优先于旧助手的固定四核/十二核门槛；科学方法、基组、电荷、多重度、网格、
VeryTight、原案例用途、两个稳定性修复循环及“采集不等于科学通过”的规则不变。

**上线边界：这里的 Python 模拟测试与 Windows Job Object 测试不等于真实 Gaussian 重启验收。
在原机完成小案例原生中断/恢复验收后再迁移昂贵作业。不要热覆盖正在执行的旧协调器，
不要擅自杀掉当前 OLD-018，也不要把旧菜单的 6 当作本版本的 6。**

## 解决的问题

原来的 `initial = Opt + Freq` 是一个不可拆的暂停单元；v2 分成 `opt-00`、`freq-00`、
`stability-00`。需要修复波函数时再接 `opt-01/freq-01/stability-01`，最多到 02。
原子仍是 SP；扫描/固定几何仍是 Force，不被悄悄优化。

中断方式有明确区别：

| 操作 | 已运行作业 | 后续行为 |
|---|---|---|
| 4 / `pause` | 当前子阶段继续到正常结束并转换/校验 | 暂停派发；协调器保留 |
| 6 / `resume` | 不打断活跃作业 | 立即补足可用槽位，不再丢失退出的 worker |
| 8 / `interrupt` | 终止本协调器拥有的 Windows Job 整个进程树，等待退出 | 冷复制/校验 chk，记录 interrupted，保留协调器 |
| 9 / `shutdown` | 当前子阶段收尾 | 协调器退出，供迁移/升级/重新选择队列 |
| 0 | 仅菜单退出 | 不隐式停止计算 |

立即中断不是 Gaussian 内部的事务式保存。当前尚未写入检查点的 SCF 工作可能丢失，
断电/硬杀时检查点也可能损坏。所有恢复先做 `formchk`、格式和原子/电荷/自旋/电子数/
基组维数身份检查；不以“存在 job.chk”充当恢复保证，坏文件保留并报错。

### 恢复策略

- **优化未完成**：从通过校验的冷检查点生成 `Opt=(Restart,VeryTight,MaxCycles=512)`，
  读取已保存的优化状态与模型化学；换核数不要求重新从最初的分子几何开始。
  若 Gaussian 拒绝内部优化历史，进入 `needs_review`，不静默退回初始结构。
- **优化已完成、解析频率中断**：复用已完成优化的检查点，仅重做解析频率。
  v2 不伪称解析 Hessian 可由 `.chk` 从任意百分比恢复，不生成错误的 `Freq=Restart`。
  `# Restart` 的解析频率续算需要保存且一致的 RWF；本版不自动使用该高级路径。
- **稳定性/SP/Force 中断**：以有效检查点的几何、实际 R/U 方法和轨道作 Guess=Read 重新执行该子阶段；
  不宣称从稳定性本征求解器某个 iteration 无损恢复。
- **Gaussian 已正常结束、协调器在转换/记账前退出**：校验并收集已完成结果，不再重跑它。
- **整个候选已经采集**：跳过，但仍保留科学裁决、竞争电子态和环境验证待办。

每次尝试新建 `attempt-NNNN`，保存输入、日志、raw chk、校验后的 chk/FCHK、资源分配、
起点文件哈希、恢复策略、真实进程树 CPU/墙钟时间。旧 attempts 不覆盖，不自动删除 RWF/CHK。
默认单次最长 72 小时；超时中断并暂停，不再采用“24 小时到期就把整个案例不可恢复地判死”的行为。
可将 `timeout_hours` 设为 0–168 范围内的正数；修改不会自动热重配正在运行的 Gaussian。

## 核心与内存

从 Windows 物理核心拓扑读取数量，不把 SMT 线程当物理核心。

| 物理核 | 预留给系统的预算 | 所有受管 Gaussian 合计上限 |
|---:|---:|---:|
| 16 或更多 | 至少 2 | 14 |
| 12 | 2 | 10 |
| 8（用户 98x3 的典型情形） | 2 | 6 |
| 6 | 2 | 4 |
| 4 | 2 | 2 |
| 2 / 1 | 1 / 0 | 1 |

少于 16 核标记为 `98x3-small-host-policy`；这是用户的机器档位名称，不是假装识别了 CPU 型号。
多处理器组机器暂时明确拒绝，不能把第一组当成整个多路服务器。

资源需求从本次 **REF-001 冻结 BSE 基组**的收缩壳层推导 NBasis，不沿用原库小基组数字。
初始参考建议核数：<=64/160/320/480/640/900/1200/>1200 个基函数分别为
1/2/4/6/8/10/12/14；大于 160 且开壳层再加 2，最高 14。
两个待发任务做整数公平分配，总和不超预算；单独的大任务（>=400 basis）可借满空闲预算。
例如 16 核机两个大任务可以 7+7，大任务与小任务可以 10+4；8 核机对应总和最多 6。
这是一套保守启发式，不是已经测过的最优加速比，不保证 14 核比 4 核快 3.5 倍。

**核数在启动一个子阶段/恢复 attempt 时选择，不会把正在跑的 Gaussian 原地从 4 核变成 14 核。**
已经在跑的长阶段不会仅因旁边空出核心而被强杀重启。

保留每作业 `%mem=24GB`、进程树硬上限 32 GiB；最多两个并发。
总预算同时受 128 GiB 和可见物理内存减 16 GiB 约束，启动前另检查可用内存余量。
因此至少需要 48 GiB 可见内存；它不会因为 8 核就错误申请两份超过本机容量的内存。
磁盘空余低于 `disk_reserve_gib`（默认 20 GiB）时暂停新派发；不自动删除证据腾空间。

### Gaussian 16W A.03 的 affinity 例外

冻结进程实现已记载：该 PGI 运行时继承缩窄 affinity 后会在读输入前崩溃。
因此本版本维持完整启动 affinity，以 `%nprocshared`、线程环境变量和 Windows Job CPU-rate
硬上限限制整个后代树。预算从物理核推导，但**不是“独占/硬绑定某 14 个物理核心”**。
Windows 可把线程调度到不同 SMT sibling；不要把预算、线程数和物理绑定混为一谈。
同机 v2 全局锁 + 旧版同路径 `jobs/supervisor.lock` 防止多个协调器叠加运行；
启动前发现其他 Gaussian/link 进程则拒绝。预算不能管住随后由不合作的外部程序另开的计算。

## 首次使用与旧结果迁移

在完整仓库 checkout 中运行，保持 `runtime_v2` 与 `paused-20260906` 是同级目录。
不要只复制 resume.py 一个文件，也不要覆盖旧 Resume/tools/runner-source。

1. 复制 `config.example.json` 为 `config.json`，核对已有 `data_root`、`work_root`、Gaussian 路径。
   示例正是当前 `F:\CalChem\COV\Resume\data`、`...\work`、`F:\CalChem\Gaussian\Gaussian 16 W`。
2. 双击 `START.cmd`，先选 **1**。它验证 273 个输入、冻结源码、Gaussian 二进制身份及物理核预算；不启动计算。
3. 先完成小案例真实 Gaussian 验收。旧 coordinator 和 Gaussian 全部退出后，选 **2**。
   导入器取得与旧版相同的 OS 文件锁，冷复制旧结果/检查点并重验 chk/FCHK 的一致性。
   导入到 `work/runtime-v2/`；原 `work/jobs/REF-001`、归档 data、暂停证据不修改。
4. 再选 **3**。默认选择全部未完成候选；也可限定数量。资源自动分配，不再手填 workers=2/cores=4。
5. 以后 **4 暂停 → 6 恢复** 可以在同一个尚存活的协调器中操作，空闲 worker 会立即补回。
   已经选 9 或协调器已退出时，6 只改变控制状态，仍需 3 重新启动协调器。

直接命令（在此目录，Python 3.12）：

```powershell
py -3.12 -X utf8 resume.py check
py -3.12 -X utf8 resume.py import-legacy
py -3.12 -X utf8 resume.py run --case OLD-018 --case OLD-019
py -3.12 -X utf8 resume.py run --limit 273
py -3.12 -X utf8 resume.py pause
py -3.12 -X utf8 resume.py resume
py -3.12 -X utf8 resume.py interrupt
py -3.12 -X utf8 resume.py status
py -3.12 -X utf8 resume.py shutdown
```

审查后重新尝试显式失败案例（协调器必须退出），必须给出理由；旧错误记录归档，命令本身不启动计算：

```powershell
py -3.12 -X utf8 resume.py retry-reviewed --case OLD-018 --reason "已审查日志并修复外部磁盘空间问题"
```

这不是跳过科学失败的捷径。重启历史被 Gaussian 拒绝或 chk 损坏时，需要原生检查/另一个经审查的起点，
不能通过删文件、改标志或者无限 retry 冒充成功。

## 状态与估时

`status.json` 有协调器 heartbeat、活跃阶段、分配核心、进程树 CPU/墙钟时间；`status` 输出
最新三组完整 Maximum/RMS Force/Displacement、阈值与通过情况，以及真实 NBasis、优化 step。
`running` 的旧文件值不是进程仍在运行的证明；文件变大也不是科学收敛证明。
不再输出“SCF 95%”“能量下降就快到极小值”或把 `Step / MaxCycles` 当完成百分比。
阈值倍数按日志打印的舍入阈值计算，仅用于诊断，不转成预计剩余小时数。

## 回归与原生验收

```powershell
py -3.12 -X utf8 test_runtime.py
```

覆盖物理核/SMT区分、1–14预算组合、8核主机、分阶段输入、不改变科学条件、R→U、两次修复上限、
坏/缺失检查点、资源变更后 Opt=Restart、解析频率重算但优化复用、旧数据不变与重复导入、
撤销暂停填回 worker、受控失败不无限重试、最新完整收敛表、全部 273 个真实冻结输入身份/基组解析。
Windows CI 另外运行真实 Python 子进程树的句柄/退出/中断测试，不使用 Gaussian 授权程序。
CI 不能替代原机上的 Gaussian 原生恢复测试，尚未有这些测试结果时必须明确 pending。

Work 原生验收应在独立临时目录完成：小优化中断→校验→Opt=Restart；不同核数恢复；
优化结束后频率中断→只重算频率；Stable 阶段重放；双作业 pause/resume/interrupt；
进程树无遗留；无误杀外部进程；单核/8核/16核预算；完整与恢复结果的能量、几何、频率、S²对照。
不要为了跑测试放松 VeryTight、缩小正式基组或覆盖正式科学证据。

## 技术依据（2026-09-08 核对）

Gaussian 官方站点本次返回 502；重启语义对照其手册的官方支持商镜像：
- https://conflex.co.jp/gaussian_support/opt.php （Opt Restart 读取检查点模型化学与优化历史）
- https://www.conflex.co.jp/gaussian_support/freq.php （数值 Freq Restart 与解析频率 RWF Restart 区分）
- https://learn.microsoft.com/en-us/windows/win32/api/sysinfoapi/nf-sysinfoapi-getlogicalprocessorinformation
- https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-jobobject_cpu_rate_control_information
- https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-updateprocthreadattribute

依赖的原验证定义：`../paused-20260906/runner-source/gaussian_reference_inputs.py`、
`gaussian_reference_batch.py`、`gaussian_rebuild.py`、`validation_process.py`、`identity.json`。
