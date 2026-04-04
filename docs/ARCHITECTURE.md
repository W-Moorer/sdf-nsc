# 架构说明 (ARCHITECTURE)

## 当前架构与隔离决策 (第二阶段)
在本项目中，我们选择将 Project Chrono 当作独立的第三方库通过 CMake 的 `find_package` 进行挂载（见 `CMakeLists.txt`），安装于 `_deps/chrono-install` 下隔离。不直接入侵 Chrono 源码。

同时我们在代码实现中引入了 **Backend 抽象层** 概念：
- **`include/platform/backend/IRigidSystem.h`**: 底层不可见具体依赖的纯虚概念，用于对上层暴露物理刷新和获取状态能力。
- **`src/backend/chrono/ChronoRigidSystemNSC.cpp`**: 专门用于封装并代理 `chrono::ChSystemNSC` 对象的组件。此时 Chrono 的对象生命周期、接触材料等底层概念完全属于平台 Backend 的内部私有状态。
- **`src/main.cpp`**: 最薄的胶水层，负责启动 Backend 并在顶层控制外部循环和基础数据输出。

> 说明：原脚本设计时输出目录期望是 `_build/platform`，由于在 Bash 子系统调用等原因，现在的实际主工程输出路径位于 `_build/project` 下。