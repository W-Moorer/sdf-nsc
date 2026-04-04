# 平台演进 Roadmap

## Phase 1: Engine Initialization (Completed)
- 探索 Chrono C++ 环境架构
- Bash 脚本自动拉取源码配置
- CMakeLists 的组织：生成 `minimal_chrono_nsc`

## Phase 2: Baseline Extraction & Refactor (Completed)
- [x] 将硬编码的魔法几何和重力约束剥离 (Target: `baseline_drop_nsc`)
- [x] `ChCollisionSystemBullet` 验证：修复 Y轴连续穿模下坠问题。
- [x] 抽象化 `DropCaseConfig` 参数数据与模型控制。
- [x] 增加命令行 (`main.cpp` 支持 `--dt`, `--T`, `--ball-height`)

## Phase 2.5: Reference Comparison Framework (Completed)
- [x] 设计支持可选列的 Reference csv 标准契约。
- [x] 编写轻量 C++ Validation 接口 (`ReferenceCurve`, `ErrorMetrics`) 解析数据及线性插值时间轴。
- [x] 指定 `BaselineDropCase` 的明确输出与 Console 总结标准。
- [x] Python 错误评估绘制脚本 `compare_drop_case.py` 增加。
- *(注：当前还没有真正导入 RecurDyn 数据，但接口已准备好。后续只要把商业软件导出的 CSV 放到 `data/reference/`，就可以直接通过 python 进行曲线对接验证与绘图。本迭代只建立对比框架，不引入 SDF / 大规模的 Block Friction。)*

## Phase 3: SMC Expansion & Geometry (Upcoming)
- 支持更精细的网格导入并对比 NSC 与 SMC 接触材质的系统表现。
- 引入 Block Friction 测试，即实现 `baseline_block_friction` （规划中，暂未执行）。
- 开始对接外接 SDF 形变渲染或深度学习碰撞补偿接口。
