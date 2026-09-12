# Work 持久交接：REF-001 快速暂停版 2.1

维护 O1dDing/CUDA-Orbital-Visualisation 的 test/fchk-validation-native。
先检查远端最新 HEAD、本地未提交改动和实际运行进程，读取 validation/README.md、
runtime_v2/README.md、FAST_PAUSE.md、resume.py、windows_job.py、windows_pause.py、
fast_checkpoint.py、native_acceptance.py、policy.py 及两份测试。
将本文件的规则写入 Work 的项目持久交接记录，不只在聊天中说“记住了”。

此次用户要求解决的是点安全暂停后等待数小时。新方案必须区分：
1. 4/hold 快速内存暂停，6/resume原进程原地继续。不等待整个优化，内存仍占用，
   不可关机、重启或关闭协调器，不能把该状态称为已落盘保存。
2. 10/pause 等通过本机验收的优化检查点段或当前子阶段结束后停止派发。
3. 8/interrupt需SAVE确认，只终止所属Windows Job树，等待全部后代退出后冷复制/刷盘/校验。
   解析Freq保留命名RWF/INT/D2E/CHK；Opt优先保存较小CHK及输入日志，原RWF不删。
   保存回执和formchk可读不是恢复历史有效性的绝对保证；未写入工作可能损失。
4. 9/shutdown收尾后退出；如果当前RAM暂停，会先恢复去完成边界，不是立即关窗。

快速暂停只操纵明确属于本Job的线程真实句柄，恢复只撤销自身增加的suspend count，
不得按进程名taskkill、根据陈旧PID误杀外部作业或使用未审查的Nt接口。
保留完整启动affinity，避免已记录的Gaussian16W A.03 PGI初始化缺陷。
72小时按active elapsed计时，暂停时间剔除；期限到达默认RAM hold，不再砍掉长任务。

Opt步级磁盘保存采用事先声明的%KJob L103 2，再冷校验与Opt=Restart，绝不能把
检测日志行后抢时间强杀当成安全边界。KJob格式、实际优化进展和结果一致必须通过原机验收；
只有native-capabilities.json中匹配当前runtime、Gaussian全套link/DLL指纹、方法和证据哈希的能力才能开启。
解析频率从不可变冷快照复制到新attempt，使用#p Restart；不使用数值频率的Freq=Restart。
RWF不可用/未通过验收时needs_review，不静默重算整个Freq。显式审查后选择replay才允许频率重放。

新入口仍在validation/runtime_v2/START.cmd；运行证据独立到work/runtime-v2-fastpause。
旧v1/v2工作目录和冻结源码、输入、身份、归档全部保留，不热覆盖旧进程。
OLD-018超时及OLD-019已完成initial等历史资料，以validation/ref001-progress-20260909为定位线索，
但必须实读本机并核验。旧timeout的有效salvage仍需明确reviewed导入，不得仅删review后fresh start。
尤其不得用OLD-018的唯一CHK做首个破坏性测试；原始SALVAGE保留只读。

资源：按物理拓扑，不按SMT。16+物理核总预算14；12核10；8核6；6核4；4核2；2/1核1。
<16核按98x3小主机档位称呼，但以真实检测值分配。单任务1–14，最多两个并发，总和不超预算。
真实本轮REF基组NBasis、开壳层和空闲资源决定分配，核数只在新attempt/段启动时改变。
每任务%mem24GB、Job树32GiB，保留系统内存和磁盘余量，协调器全局锁与旧目录锁同时使用。
模型、基组、电荷/自旋、VeryTight、SuperFineGrid、固定几何用途、两次稳定性修复上限不变。
candidate_collected不是scientific_pass；SCF小误差、日志变大都不是完成百分比。

执行：
先只读检查、备份并运行test_runtime.py和test_fast_pause.py、原冻结恢复测试。
再在本机、无其他Gaussian且取得相同OS锁后运行：
F:\Dev\Python312\python.exe -X utf8 resume.py native-acceptance
它在独立临时目录自动对水和NO的副本做baseline、RAM暂停继续、KJob优化分段/换核数、
解析Freq中断/RWF Restart并对比能量、原子间距离和频率，产生真实日志及能力回执。
某项未实际触发或失败不算通过；不能编辑capability=true绕过检测。
CI的Windows Python父子进程测试和模拟测试不等于Gaussian原生通过。

如果原机测试暴露问题，保留失败日志、修复并加回归后重新原生验收；不得放松科学阈值换通过。
验收通过后再部署正式入口、做reviewed迁移；不擅自强杀正在跑的旧任务、不删除锁、
不重新算已有效完成的阶段、不擅自发COV release或改main/feat分支。
将最终commit、实际入口、配置/数据路径、能力回执、测试结果、迁移结果和剩余限制持久记录。
