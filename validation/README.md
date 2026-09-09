# 验证项目入口

## 2026-09-09 COMMON-ROOT-008 新全量验证

[COMMON-ROOT-008 稳定候选检查点](common-root-008-20260909/README.md)已保存：
当前源码8e206e1的完整专项控制、固定40预检和流程门检查通过，新的273例采集已启动。
原始轮次和复核队列保持冻结；中断后先核实实际终态，只接续缺失阶段。
该记录不表示全量复核或完整科学验收通过，未启动Gaussian。

## 2026-09-09 当前接续检查点

[REF-001 增量归档与恢复记录](ref001-progress-20260909/README.md) 已保存并上传，
[远端校验回执](ref001-progress-20260909/upload-verification.json) 确认 4 个新增资产，
合计 959,732,968 bytes，大小及 SHA-256 全部一致；旧数据归档和本地原件保留。
现有证据确认 17 个已收集候选、37 个已完成阶段、OLD-018 超时断点及 OLD-019 已完成 initial。
原打包 manifest 的 pending 是保留的历史状态，以后续上传回执为完成依据。

最新用户顺序是归档和恢复标记完成后返回 COV COMMON-ROOT-008，
不要求先立即重启 Gaussian 或完成 runtime v2 迁移。Runtime v2、物理核资源政策和科学精度要求仍有效，
真实 Gaussian 原生重启验收仍待完成。详情及后续验证门槛见 [handoff.json](ref001-progress-20260909/handoff.json)。


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
