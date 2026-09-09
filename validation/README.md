# 验证项目入口

## Gaussian 续算执行层

新的执行入口是 [REF-001 Runtime v2](runtime_v2/README.md)，位于
`validation/runtime_v2/START.cmd`。资源预算和中断恢复采用该目录的新规则：
按实测物理核心分配，16 核及以上最多合计 14 核，8 核通常最多 6 核；
Opt/Freq/Stable 分开，已校验检查点恢复、暂停后恢复空闲槽位。
该执行层上线前仍须完成真实 Gaussian 16W 的原机中断/恢复验收，不能把模拟/CI 通过当成原生通过。
给本地 Work 的完整交接指令见 [WORK_PROMPT.md](runtime_v2/WORK_PROMPT.md)。

`paused-20260906/` 是不可改写的旧版冻结快照。保留它的输入、源码、身份及证据；
不要在原目录直接替换 runner 或删除锁文件。当前正在运行的旧版进程不会自动获得 v2 功能。
新证据记录在本地 `work/runtime-v2`，旧 `work/jobs/REF-001` 只读导入。

## 科学规则与其他验证

科学要求仍参考 `../docs/FCHK_VALIDATION_AGENT.md`、`../docs/native-validation.md`，
以及各日期目录内的具体冻结契约、真实证据和待办。此次新用户资源授权取代旧固定 4 核、
12 核门槛等执行限制，但不改变模型化学、基组、电子态、网格、收敛精度和科学验收标准。

`candidate_collected`、程序正常退出、界面能够显示、自动测试通过，都不等于完整科学验收通过。
