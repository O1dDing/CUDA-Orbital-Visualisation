# REF-001 增量归档与恢复检查点 — 2026-09-09

本包保存计算证据，不是软件发布或科学验收通过声明。扫描时没有 Gaussian/link 进程；已实际取得并释放旧工作目录的 OS 互斥锁，没有删除锁文件或撤销 PAUSE。

| 核验结果 | 数量/状态 |
|---|---|
| 已采集参考候选 | 17：OLD-001 至 OLD-017 |
| 已完成且重新核验的阶段 | 37 |
| OLD-018 | initial 超时，可恢复候选；优化未完成，频率未开始 |
| OLD-019 | initial 的优化和频率均已完成，阶段间暂停，下一步是 stability |
| 尚未启动案例 | OLD-020 至 OLD-273，共 254 |
| 本次启动计算 | 0 |
| 完整科学验收 | 尚未完成；candidate_collected 不等于 scientific_pass |

OLD-019 的 result.json 仍写 running，但已完成阶段、Normal termination、FCHK、checkpoint、实际输入身份和增量调度器的 paused_between_stages 记录互相支持。保留原状态原字节，以 case-status.json 给出本次核验后的解释。37 阶段均核对了输入哈希、FCHK 哈希、真实参考身份、checkpoint 存在及当前哈希、Gaussian 正常结束和 formchk 成功记录；没有为重复验证而重新计算。

## 内容和去重

manifest.json 列出完整工作目录和 salvage 文件的来源、时间戳、大小、SHA256及恢复定位。新增的102份独立内容打成3个 ZIP，包含 OLD-013 的后续稳定性修复链、OLD-015/016/017、OLD-018 超时现场与 RWF/CHK/log、OLD-019 已完成 initial，以及新的状态和用户原始要求。

与旧归档字节一致的文件记录 existing_archive，引用 [旧续算资料资产](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906) 内的具体 ZIP、内部路径及资产 SHA256；本次已通过 GitHub 返回的资产 digest 核验旧 manifest 和全部23个分卷的身份，不重复上传。新文件之间的同内容副本记录 same_content_as，但保留每个原始路径。

本地原件、原 data、旧归档和冷备均保留。Python 缓存和临时 OS 锁不作为可恢复科学数据；不要从归档恢复正在持有的锁。许可证软件本体未打包。

## OLD-018 的恢复边界

原始 job.chk 与 salvage/job.chk 的 SHA256 均为 `12c50a01ca0d2b37b1f1234df94ea7d0b2c62553ef605a1108212238de08087e`，27,983,872 bytes。已有 salvage FCHK 重新通过冻结参考的原子、电荷、多重度、电子数、方法和基组身份校验：UPBE1PBE，NBasis=502。历史交接记录过成功 formchk；本次复用了该 FCHK，没有重复转换。恢复前仍须冷校验并在小案例完成真实 Gaussian 16W 中断/重启验收。

超时是旧运行器86400秒阶段期限，记录exit_code=124、timed_out=true、SCF failure=false。不能把超时称为已优化完成。以后优先从这个有效断点用 Opt=Restart 建立新 attempt；若Gaussian拒绝优化历史，保留拒绝日志进入needs_review，不静默从原始结构重算。

MaxCycles调查：实际 job.gjf 和日志回显都是512，Gaussian展开的路由明确出现 `1/6=512`，因此不能归因于生成器没写入关键字。Berny同时显示maximum=150；二者不一致的内部原因仍待专门原生检查，本次不猜测、不改变科学参数。原日志完整保留。当前退出发生在step46的运行器时限，证据不支持声称它触及150步终止。

## 最新有效资源和阶段规则

Runtime v2 与本地参考缓存改动已保留并整合。当前主机16物理核/32线程，Gaussian合计最多14核、最多2个作业；核数按真实REF基组规模、开壳层及剩余预算在新阶段/attempt时分配，不热改运行中作业。每任务输入24GB、Job树32GiB，并核验总内存、可用内存、磁盘和同机互斥。模型、基组、VeryTight、SuperFineGrid、稳定性及适用的两次修复上限保持。

Opt/Freq/Stable分开保存。优化中断优先校验后的Opt=Restart；解析频率中断保留优化，只重做频率；Stable/SP/Force从有效检查点保存的几何、实际R/U方法及Guess继续该阶段。不能把chk说成可恢复任意解析Hessian进度。旧v1菜单没有v2的新语义。

本机runtime_v2自动测试已成功，含真实Windows Python进程树控制；旧只读check确认273输入和37阶段。真实Gaussian重启验收仍为pending，不能据此把v2称为已完成原生验收。本次归档不部署、不启动或迁移计算队列。

## 恢复与接续顺序

先保留并核验现有本地文件。需要异机恢复时，在新的隔离目录按manifest定位各文件：新ZIP取archive字段，旧内容取existing_archive的ZIP与内部路径，同内容路径取same_content_as；逐项核对bytes与SHA256。先还原旧23卷的data与新工作快照，绝不覆盖唯一现有工作目录。原始证据里的绝对路径保持原样；实际迁移只重定位新工作副本，重新执行冻结输入/Gaussian/checkpoint身份门槛。v2只能在小案例原生验收及旧协调器退出、取得同一路径锁后冷导入，原work/jobs/REF-001不改写。

用户最新澄清的顺序是：全部归档上传并做好恢复标记后，立即回到COV主线；不以立即续算OLD-018/019或完整迁移验收作为返回前提。“撤回runtime v2或回退计算修改”是已纠正误解，不执行。

COV从COMMON-ROOT-008合并草稿接续。NUM-FIX-007C的273例已完成取证及限定范围裁决，核验后复用。81项准备控制早于六种单位标签修复；先验证新增变化、真实调用和渲染覆盖。局部投影/标签作用域、骨架螺环、能隙/配体先验和视区比例集中完成，专门控制通过后跑固定40例预检，经workflow_gate和实际回执再冻结正式273。仅检查器改变且旧原始证据完整时另建review复用，不重新取证。独立缓存只复用身份完全匹配的参考，不缓存COV结果或通过结论。最终仍须原273加至少50个新外部分子的适用完整验收。

归档资产上传核验见upload-verification.json；本地完成不等于远端上传成功。
