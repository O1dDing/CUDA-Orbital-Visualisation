# Gaussian adaptive physical-core scheduling — 2026-09-13

本次交付进入 `feat/fchk-primary-interface`，在 PR #3 已有主程序代码上增加统一的 REF-001 自适应调度。只选择性引入当前实际使用的运行时和 Helper 依赖，没有合并整个验证分支。原始输入、计算产物、CHK/FCHK/RWF、暂停记录与生产者身份保持原位。此次没有启动正式273例批次或 external50。

## 273例完整资源表

- [全部273例 CSV](all-273.csv)
- [全部273例 JSON，含每个实际输入及NBasis日志的路径和SHA](all-273.json)
- [全部273例可读表及完整重任务表](all-273.md)
- [全部31个大于5核的重任务 CSV](heavy-over-5.csv)
- [汇总、历史实际核数与原生验证结果](audit-summary.json)

| 主机 | 总预算 | 小档 | 中档 | 重档 | 同时任务 |
|---|---:|---:|---:|---:|---:|
| 16个物理核 | 14 | 4核，203例 | 5核，39例 | 7核，31例 | 最多2 |
| 8个物理核 | 7 | 2核，203例 | 3核，39例 | 4核，31例 | 最多2 |

阈值统一为 NBasis <300、300–449、≥450；8核档位为16核档位除以2后向上取整。两组计数各自合计273。所有已完成案例也参与分类，但分类不会重新调度它们。名字取自冻结输入标题，分子式按原子序数计算；标题不是重新进行化学身份鉴定。

109例具有同一REF方法的Gaussian日志NBasis，均与冻结基组计数相符。其余164例严格读取 `basis.gbs`，逐原子累加收缩轨道壳的维数：S=1、P=3、SP=4，纯球谐壳为2l+1，Cartesian壳为(l+1)(l+2)/2。弥散壳按实际基组计入，ECP势项不冒充轨道基函数。解析器对未支持的广义收缩格式拒绝猜测；当前273份冻结基组均通过。旧的小基组FCHK中的NBasis没有用于分类，原子数也没有充当NBasis代理。

## 为什么使用这些档位

REF-001使用含25%精确交换的PBE0与D3(BJ)，基组为def2-TZVPP或带弥散函数的def2-TZVPPD。多数案例还包括VeryTight优化、CalcFC初始Hessian、Freq解析Hessian、VeryTight/XQC SCF和SuperFineGrid。基函数维数影响矩阵、交换积分与响应计算的工作量，因此用于公共档位；重元素/ECP、弥散函数、开壳层和已有收敛失败作为逐例理由保留。这些特征不能证明更多核必然带来最优速度，也不能成为削弱科学设置的理由。本次未用运行时暂停失败推断分子的SCF收敛性。

OLD-108的冻结标签为anthracene，分子式C14H10、24原子，电荷0、自旋多重度1。实测NBasis=574，PBE0-D3(BJ)/def2-TZVPP、5D7F、VeryTight Opt和SCF、CalcFC、Freq、SuperFineGrid均已核实。按共同阈值，16核机器分配7核，8核机器分配4核；没有OLD-108专用分支。其旧实际运行输入曾使用8核，历史保持不动。

## 真实执行链与历史4+4

`validation/paused-20260906`保存的是冻结旧执行链。原始参考输入是 `%nprocshared=4`，旧 `gaussian_reference_batch.py`的线程池最多两个任务，`resume_reference.py`的每槽分配也固定为4，因此双槽实际是4+4。这与机器后来的14核总预算是不同层面的设置；预算14不会自动改写冻结4核输入。

`validation/ref001-progress-20260912`是当日归档，记录90例候选已收集、2例部分完成、2例待复核。它是历史快照，不是当前运行器的资源策略或当前总进度。历史 `work/runtime-v2-fastpause`采用1/2/4/6/8/10/12/14偏好、剩余量分摊和单个重任务放大到预算的规则。实际保存的输入中，284份分别为1核73份、2核93份、4核90份、6核28份；所以过去确实存在大于5核运行，不能概括成始终4核。

本次迁移前的实际入口是 `Resume/runtime/active.json` 指向的统一2.2代码包，Helper、菜单和Resume共用它。该运行时在8物理核时得到6核预算，且存在单个重任务放大到14、按剩余量降档的问题。当前统一目录保存的10份实际运行输入中，4核5份、6核2份、8核3份。所有逐输入出处与哈希在JSON表中，数字不是根据界面推测。

本次把 `policy.py → Engine.run → calculate_reserved → Gaussian/formchk run_tree` 接到一份持久资源账本。动态 `%nprocshared` 仅写入独立attempt输入；冻结资源头仍为4。已有计算记录的生产者身份、输入哈希和科学路径不会被新策略改写。

## 准入、归还、暂停与崩溃

申请完整档位必须同时满足“已分配总数+本次申请≤预算”和“并发任务数<2”，然后才提交工作线程。16核允许4+4、4+5、5+5、7+4、7+5、7+7；8核允许4+3、拒绝4+4。队首重任务等不到完整档位时，后来的小任务不能绕过它；已完成阶段转到等待队列末尾。

一份阶段预留覆盖Gaussian整棵进程树及其串行formchk等辅助工具；独立原生探针也必须预留。树的子租约禁止并行复用同一个阶段名额。RAM暂停保留原核数和任务名额，记录保留原因；恢复同一树不会把资源释放后悄悄超卖。正常退出、失败、取消和阶段暂停只有在整棵树退出后才能归还。

账本记录申请/分配/预算、并发任务数、策略版本、参考和运行时身份、输入所在目录、PID及创建时间、Job名称与分配/树启动/确认退出/归还/暂停事件。新子进程通过Windows `PROC_THREAD_ATTRIBUTE_JOB_LIST` 在创建时直接纳入Job，关闭了创建后再纳管的崩溃窗口。重启核对PID创建时间和Job存活，仍活着的树继续占预算；不凭心跳过期删预留。所有事务跨线程及跨进程序列化。手工在托管执行链外启动的Gaussian仍属于外部进程，入口发现后会拒绝派发，不会擅自接管或终止。

## 验证与实际边界

回归覆盖全部16核组合、8核阈值与4+4拒绝、第三任务拒绝、线程/进程竞争、原生创建时崩溃、PID复用、存活孤儿、失败清理、内存暂停/恢复和队列顺序。还保留既有科学身份、冷保存、暂停回滚和Helper不自动计算的检查；最终测试记录见 `tests-summary.json` 与 `tests.log`。

隔离的水/NO在4核档位执行8个真实Opt/Freq阶段：基线与暂停恢复后的能量、几何和频率一致，两种R/U方法均通过该新运行时身份的RAM暂停验收。详见[native-pause-report.json](native-pause-report.json)。L103优化分段与解析RWF恢复继续保持未验收门槛。

OLD-108使用原冻结科学设置的独立副本，以7核运行60.11秒后按计划停止；NBasis=574，Gaussian输出确认7线程，OMP/NCPUS/OMP_THREAD_LIMIT/MKL均为7，整个Job共计19个进程，累计385.875 CPU秒，平均6.419个CPU时间当量。Windows Job CPU rate为2187/10000（在32逻辑处理器机器上约7个CPU时间当量），所有后代确认退出。输入、日志、环境和采样见[native-evidence](native-evidence)，冷保存CHK/RWF及完整原始回执留在本机隔离检查目录。该smoke不是已完成优化，也不是性能最优或连续满7核的证据。

物理核数取自Windows硬件拓扑，16个物理核没有被32个SMT线程替代。预算和档位按物理核决定，运行时叠加Gaussian线程设置与整个Job的CPU时间硬上限；为保持Gaussian16W A.03兼容性，完整启动亲和性保留，未声称把进程固定或独占在7个物理核上。8核场景做了策略和竞争测试，本机没有8核硬件实测。参考：[Windows Job创建属性](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-updateprocthreadattribute)、[Windows CPU Sets的软亲和性说明](https://learn.microsoft.com/en-us/windows/win32/procthread/cpu-sets)。

## 部署与保留

新消费身份为 `d9f48dc70e1edab1a9150bb55eb0d661e4c33f6a7996edfceb086ca711f95db4`。迁入 `Resume/work/runtime-v3-adaptive-d9f48dc70e1e` 的345个收集阶段、105例已完成候选、2个恢复来源和2个待复核标记均核对过旧记录及原始文件SHA；迁移不启动计算，不复制旧native能力回执。详见[migration.json](migration.json)。

安装使用新的不可变 `runtime/releases/<bundle>` 与身份迁移，最后只切换共享配置和active指针。原2.1/2.2目录、冻结归档、科学输入、结果与OLD-018-SALVAGE保留，旧代码包仍可用于按回执恢复。打开新Helper默认只读、手动控制关闭，不自动恢复或开始队列。最终部署回执位于 `Resume/runtime/deployments/adaptive-physical-20260913.json`。PR #3保持未合并，PRE11和其他Release没有改动。
