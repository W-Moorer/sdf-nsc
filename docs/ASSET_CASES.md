# Asset Cases 清单与标准

本文档汇总了当前工程中 `assets/` 目录下扫描到的所有物理案例 (Cases)，并明确保留原有的 Chrono 最小验证基线。

## 扫描到的 Asset Cases 列表

### 1. `cam`
- **状态 (Status):** `ready`
- **描述:** 凸轮机构案例
- **建模文件 (.rmd):**
  - `simple_cam.rmd`
- **几何文件 (.obj):**
  - `models/cam_body1.obj`
  - `models/cam_body2.obj`
- **结果文件:**
  - `data/cam_data.csv`

### 2. `simple_gear`
- **状态 (Status):** `ready`
- **描述:** 简单齿轮啮合系统
- **建模文件 (.rmd):**
  - `simple_gear.rmd`
- **几何文件 (.obj):**
  - `model/GEAR21.obj`
  - `model/GEAR22.obj`
- **结果文件:**
  - `data/Gear22.csv`
- **其他文件:** 包含 Python 提取脚本与参数设定（如 `extract_objs.py`, `locked_parameters.json`）。

---

## 物理比较状态释义
系统会对所有提取的案例输出一个 JSON Manifest （位置：`data/reference/assets_case_manifest.json`）。其中记录的 `status` 表示：
- **`ready`**: 案件具有最小闭环要求的 1 个 `.rmd` 建模源文件、1 个以上的 `.obj` 几何文件、1 个以上的标定结果文件。
- **`missing_geometry`**: 缺乏 `.obj`，不可直接导入 Chrono 网格。
- **`missing_result`**: 缺乏 CSV/TXT 等比对结果，无法进行 ErrorMetrics 对比评估。
- **`unknown_status`**: 其他缺失问题。

---

## Chrono NSC baseline examples to preserve

在对接并开发商业模型（如上述的 `assets/`）之前，工程中**必须保留**并维护好的基本物理基线如下：

1. **minimal NSC sphere-ground** 
   *这是整个项目早期探索 Chrono 核心引擎、配置 CMake 所采用的核心验证模板，确保基于 `ChCollisionSystemBullet` 的底层能够不穿模运行。后续作为轻量检查用。*

2. **baseline_drop_nsc**
   *第一个被正式标准化的基线。其解耦了 `DropCaseConfig` 和业务层，自带完整的 `.csv` 输出能力，它也是对比框架 (Reference Curve / RMSE Computation) 能跑通的基准要求。绝对不允许将其删除。*

> **注意：**
> 如果后续在研究 Chrono 源码（如 `Chrono/doxygen/documentation/` 或者相关教程）中发现了适合做 NSC（Non-Smooth Contact）的示范例子（比如摩擦、堆叠等），应该以“参考映射”的方式克隆其逻辑到我们的平台上作为新的 `baseline_xxx_nsc`。切勿用其他模型直接覆盖替代它们。
