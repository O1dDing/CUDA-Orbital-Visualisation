# 从源码编译 Chemical Orbital Visualiser

[English](BUILD.md)

普通 Windows 下载包不需要 CUDA Toolkit、CMake 或 C++ 编译器。现有包适用于 RTX 50 系列；其他 NVIDIA GPU 架构可从源码构建。

## 环境要求

- 支持 CUDA 的 NVIDIA GPU 和兼容的驱动程序
- CUDA Toolkit 12.8 或更新版本
- CMake 3.28 或更新版本
- 支持 C++20 的编译器
- OpenGL 2.1 兼容上下文
- Git；CMake 用它获取 GLFW 和 Dear ImGui

默认 CUDA 架构为 `sm_120`，并包含 `compute_120` PTX。其他 GPU 需将 `CMAKE_CUDA_ARCHITECTURES` 设为目标架构。

## 选择版本

`v0.3.0` 是稳定版，`v0.4.0-pre.1` 是 NBO 预览版。以下命令在所选版本的源码目录中运行。

## Windows / Visual Studio 2022

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --parallel
.\build\Release\cov.exe .\examples\h2.molden
```

## Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/cov ./examples/h2.molden
```

## 不使用 CUDA 的核心测试

```bash
cmake -S . -B build -DCOV_ENABLE_CUDA=OFF -DCOV_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## 输入文件

COV 读取 Gaussian FCHK/FCH 和 Molden 波函数。打开 CHK 需要单独安装 `formchk`；找不到转换器时，将 `COV_FORMCHK` 设为可执行文件的完整路径。Gaussian `.log` 或 `.out` 可补充波函数信息，不能替代 FCHK/FCH。

每个输入最多 100 个原子。支持 Cartesian 和实球谐 `s/p/d/f/g` 基函数。Molden 展开后的基函数数目须与轨道系数一致。

Gaussian/NBO 文件准备见[单次作业模板](NBO_ONE_JOB.zh-CN.md)。
