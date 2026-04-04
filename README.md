# Baseline Contact Platform

## 当前进度
目前处于 **Reference 框架对接阶段**:
- 完成了第一阶段的 Chrono 下载、最简化编译安装验证闭环。
- 完成了底层 Chrono C++ 程序的初步探测、CMake 链接验证及简单的 NSC（Non-Smooth Contact）封装。
- 引入了 `platform_backend` 给外部调用的早期雏形，将小球掉落测试提取出通用的 `baseline_drop_nsc`。
- 从 hardcode 中抽象出 `DropCaseConfig` 参数控制以及 `BaselineDropCase` 面向对象模型。
- **(New!)** 建立 Reference 对比框架：让 `baseline_drop_nsc` 已经具备 reference 对比能力。目前还没有真正导入 RecurDyn 数据，但接口已准备好。后续只要把商业软件导出的 CSV 放到 `data/reference/`，就可以直接与 `baseline_drop_nsc` 输出的基准数据进行比较，甚至自动出图。

## 目录结构
- `scripts/`: 用于自动化构建、运行比对等脚本 (如 `compare_drop_case.py`)。
- `apps/baseline_drop_nsc/`: 独立出来的测试例程入口 (main) 以及命令行参数解析。
- `src/backend/`: 平台核心代码，封装底层 Chrono 物理引擎。
- `src/models/`: `BaselineDropCase` 的实现与执行逻辑。
- `src/validation/`: 引入的 `ReferenceCurve` 与 `ErrorMetrics` 对比核心框架。
- `data/reference/`: 标准对照的 CSV 参考数据（商业软件或对照结果，如 `baseline_drop_nsc_reference.csv`）。
- `data/outputs/`: 仿真输出的结果数据及生成的对比误差图表（如 `plots/`）。

## 下一步 (使用指南)
请阅读 [docs/BUILD_AND_RUN.md](docs/BUILD_AND_RUN.md) 了解如何执行构建、运行基线仿真以及进行参考比对。