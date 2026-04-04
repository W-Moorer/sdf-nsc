# 构建与运行指南

本仓库依赖 C++17 和 CMake 3.19+ 分析构建。物理后端使用 [Project Chrono](https://projectchrono.org/)。

## 1. 自动构建过程
系统内置了 Bash 脚本，可在装有 `cmake` 与 `make` 及 GCC 编译器（或者 Clang，以及 MSYS2/MinGW）的环境（推荐直接使用 WSL/Ubuntu）下全自动构建：

```bash
# 构建基础 Chrono 依赖并直接编译当前平台
bash scripts/build_project.sh
```

*(如果尚未安装 Chrono，请先执行 `bash scripts/bootstrap_all.sh` 进行自动编译依赖。)*

## 2. 运行 Baseline Drop Case (NSC)

构建完成后，在 `_build/project/` 目录下生成独立的仿真执行程序。

### 基本指令

```bash
# 使用默认配置执行
./_build/project/baseline_drop_nsc
```

### 携带参数指令

目前平台将常量固化为了面向对象的 configuration，您可以根据需求自由定义如下参数。
注意：`baseline_drop_nsc` 已经设定了规范的输出契约，如果您的程序不需要指定自定义路径，系统默认会将输出保存到：`data/outputs/baseline_drop_nsc.csv`。

```bash
# 指定 dt, 仿真时间, 输出 CSV 及小球的初始下落高度
./_build/project/baseline_drop_nsc --dt 0.005 --T 1.5 --output data/outputs/baseline_drop_nsc.csv --ball-height 4.0
```

运行结束后：
1. 会在命令行输出每一次时间步长记录的监控状态，并在结尾输出涵盖所有必要项 (`Min Y`, `Final Y`, `Final Vel`, `Had Contact`, `First Contact Time`, `Output CSV path`) 的 Summary。
2. 在指定的输出中保存包含 `Time`, `Pos_Y`, `Vel_Y`, `Num_Contacts` 必需列的高密度 CSV 以供验证框架直接读取比较。

## 3. Reference 评估对比

当 `baseline_drop_nsc` 具备输出标准数据后，您可以对接我们的验证评估框架（支持与未来外部数据软件接入如 RecurDyn），自动计算误差指标甚至出图。

### 执行指标计算与对比

```bash
# 执行比较并出图
python3 scripts/compare_drop_case.py --sim data/outputs/baseline_drop_nsc.csv --ref data/reference/baseline_drop_nsc_reference.csv --out-plot data/outputs/plots/drop_nsc_comparison.png
```

如果没有指定路径，脚本会直接搜索默认的标准输出和参考模版。如果环境拥有 `matplotlib` 库将把绘制的统计分析结果放到 `data/outputs/plots/`。

该评估框架已准备好，支持缺失 `Num_Contacts` 的曲线导入，通过线性插值强制补齐不同的时间戳，从而实现严格比较。
