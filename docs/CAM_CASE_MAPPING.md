# CAM Case 映射文档

此文档定义并说明把 `assets/cam/` 下的凸轮机构引入到我们的物理评估框架的具体信息。这是一个标志性的步骤：`cam` 作为**外部对标 Case**，也是我们将外部带有复杂接触的动力学仿真结果正式与 Chrono 系统校对对比的第一步尝试。

## 目标与定位

1. **`cam`：外部对标及真实环境引入 Case**
   - 目的在于用真实的、来自第三方软件的网格与物理记录，对接平台当前的 ErrorMetrics 和比较系统。
   
2. **`minimal NSC sphere-ground` 和 `baseline_drop_nsc`：Chrono NSC 算法基线**
   - 作为系统底层基础能力测试用（不依赖复杂的外部输入模型，可独立验证 Chrono 引擎本身的正常运作及编译器对接是否完备）。
   - **注意**：这两类基线职责不同，彼此互补。`cam` 此类外部用例绝不可用来替代或覆盖掉 `baseline_drop_nsc` 的实现。

## 使用的资产文件 (Assets Files)

`cam` 用例主要包含并依赖于以下文件：

- **`simple_cam.rmd`**
  - **角色**：第三方软件（通常是 RecurDyn 等动力学软件）的机构建模板源文件，记录了运动副约束以及动力学环境配置。
- **`models/cam_body1.obj`**, **`models/cam_body2.obj`**
  - **角色**：物理几何体模型。分为处于驱动的驱动凸轮本体和被驱动件（如推杆）。未来需要交由 Chrono 系统的 Mesh 导入器进行渲染。
- **`data/cam_data.csv`**
  - **角色**：原版（第三方）软件产出的仿真参照（Reference）真值。我们的对比框架正是使用这份文件作为计算 RMSE 与绝对误差的 Target。

## 映射现状 (Mapping Status)

### 1. 目前已完全接入的信息：
- **真值结果对比框架 (`cam_data.csv`)**
  能够通过新建的 `compare_drop_case.py` 脚本并且携带参数 `--case cam` 和 `column_maps/cam.json`，无缝过滤解析复杂的头字段名列标并直接生成比对评估误差。
- **几何与位置配置结构** 
  利用了新建的 `case_configs/cam.json` 将这些组件位置进行系统层面的声明定义。

### 2. 尚需人工配置的信息清单（TODO 依赖未来的 RMD 解析器实现）：
当前阶段，仍然因为尚未构建专门解析 `.rmd` 的 C++ / Python 层组件，一些环境变量并未做到在 C++ 项目内自动对齐加载：
- `.obj` 体质量、各向转动惯量的自动解析。
- `simple_cam.rmd` 内所含括的光滑滑动摩擦力（Smooth friction/contact material）参数设定。
- 凸轮预设的输入转速 (`Rotational Velocity`) 与重力场限制条件。

## 为什么选择将 `cam` 选做第一个接入的真实 Case？

相较于具有间歇性、高非线性接触碰撞冲量的 `simple_gear` 齿轮模组，凸轮工作表面具有**连续贴合滑移 (Continuous Sliding)** 的摩擦特征。这在检验我们在 Chrono 的接触包络及数值积分稳定性时，更能直观反应我们的系统有没有出现严重的能量漂移。连续滑轮机构曲线能给我们的 Validation 系统出具非常优美完整的误差矩阵参考。

## Compare 框架运行状态与后续计划

### 1. 当前 Compare 是否已跑通
**已成功跑通。** 现已支持直接通过 Python 脚本加载 data/reference/case_configs/cam.json 和 data/reference/column_maps/cam.json，并能够基于 ssets/cam/data/cam_data.csv 文件生成严格的误差曲线。由于当前 baseline_drop_nsc 仅是一个带属性的小球，与凸轮的位置和速度有数量级的偏差，这个巨大的 RMSE 也反向证明了真值对照的数据被正确的解析提取和对齐。

### 2. 哪些量已经成功对齐
- Time: 时间轴成功剥离 (idx 0)。不仅剥离成功，还已支持插值到 Chrono 产生的密集步长中。
- Pos_Y: Y 轴位置成功对接 (idx 10)。
- Vel_Y: Y 轴线速度成功对接 (idx 13)。
- 可基于这些核心参数成功输出 最大绝对误差 (Max Abs Error) 和 均方根误差 (RMSE) 矩阵以及趋势绘制。

### 3. 哪些量还需要后续从 rmd 或商业软件里补充
要将 cam 打造为完全可复现的 Chrono baseline case，我们未来还需要：
1. **物理模型参数引入：** 当前缺乏对两组 .obj 具体材质关联，如推杆质量、凸轮驱动转速。需要从 .rmd 解析。
2. **接触法向分析：** 凸轮不仅仅具有 Y 轴运动，我们后续可能需要引入 Pos_X、Rot_Z (绕 Z 轴旋转)、和接触力（Contact Force) 列的数据列匹配。
3. **Mesh解析器：** 真正的 aseline_cam_nsc 需要把 cam_body1.obj 和 cam_body2.obj 通过真正的 C++ Chrono 代码引入渲染以取代现在简单的下落小球。这需要编写相关的 Wavefront OBJ 解析后端。
