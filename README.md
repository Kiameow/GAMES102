# GAMES102 — 曲线插值与逼近（Qt 作业工具）

用于 GAMES102 课程作业的交互式曲线工具：在窗口中用**鼠标左键**点出红点（数据点），
一键显示**插值**（Lagrange 多项式）或**逼近**（最小二乘多项式）结果。

## 功能

| 操作 | 说明 |
| --- | --- |
| 鼠标左键 | 添加一个红点（数据点） |
| 鼠标右键 | 删除最近的一个红点 |
| 「插值-Lagrange」 | 显示 Lagrange 插值，曲线经过全部数据点（蓝色 #1f77b4） |
| 「插值-幂基」 | 显示幂基多项式插值（橙色 #ff7f0e） |
| 「插值-Gauss基」 | 高斯基函数(RBF)插值，σ 可调（紫色 #9467bd） |
| 「逼近-幂函数最小二乘」 | 显示最小二乘多项式拟合，次数可调（绿色 #2ca02c） |
| 「逼近-岭回归」 | L2 正则最小二乘，λ 可调（青色 #17becf） |
| 「显示所有」/「隐藏全部」 | 一次显示/隐藏所有算法曲线（不同颜色叠加） |
| 「清除所有点」 | 清空画布 |

> 提示：添加/删除数据点后，已显示的曲线会自动重新计算。

## 环境要求

- Windows 10/11
- **Visual Studio 2019 或 2022**（安装时勾选「使用 C++ 的桌面开发」）
- **CMake ≥ 3.21**（https://cmake.org/download/ 或 VS 自带）
- **Qt 5.15.2 (msvc2019_64)**：`C:\Qt\5.15.2\msvc2019_64`
  - 本机的 Qt 是 MSVC2019 编译的，与 VS2022 (v143) 二进制兼容，可直接使用；
  - 若你的 Qt 装在别处，构建时用 `-QtDir` 指定即可。
- **Eigen 3.4.0**：已随仓库打包在 `third_party/eigen/`（纯头文件库，MPL2 协议，
  许可证文件随附），**无需安装**。曲线算法里的线性代数统一用它。

## 快速开始

在**项目根目录**打开 PowerShell，执行：

```powershell
# 1. 一键构建（自动探测 VS 和 Qt，Debug 版）
.\scripts\build.ps1

# 2. 运行
.\scripts\run.ps1
```

常用变体：

```powershell
.\scripts\build.ps1 -Config Release     # 发布版
.\scripts\build.ps1 -QtDir D:\Qt\5.15.2\msvc2019_64   # 指定 Qt 路径
.\scripts\build.ps1 -ConfigureOnly      # 只配置不构建（方便打开 IDE）
```

### 只用 configure.bat 配置（VS2022 + v142 工具集）

如果你习惯用 Visual Studio 或命令行直接构建，可以双击 / 运行根目录的
`configure.bat`，一键完成 CMake 配置：

```
configure.bat                     # 使用默认 Qt 路径 C:\Qt\5.15.2\msvc2019_64
configure.bat D:\Qt\5.15.2\msvc2019_64   # 指定 Qt 路径
```

它会生成 `build\msvc2022-v142\GAMES102.sln`，特点：

- **生成器**：`Visual Studio 17 2022`（x64）
- **工具集**：`v142`（MSVC2019），与本机 Qt 5.15.2 msvc2019_64 完全同代，
  避免 v143 的二进制兼容性顾虑
- 之后用 `cmake --build build\msvc2022-v142 --config Debug` 或直接打开
  `.sln` 按 F5 构建运行（运行前需把 `C:\Qt\5.15.2\msvc2019_64\bin` 加入 PATH，
  或对 exe 跑一次 `windeployqt`）

> 注意：`configure.bat` 是纯 ASCII + CRLF 换行（cmd 要求），
> 编辑时请保持这两点，否则双击会报解析错误。

构建完成后 `build\msvc2022\Debug\`（或 Release）下会有一个**自包含的**
`GAMES102.exe`（windeployqt 已把 Qt/VC 运行库复制进去），双击即可运行。

### 手动构建（不依赖脚本）

```powershell
$env:QT_PREFIX = "C:\Qt\5.15.2\msvc2019_64"          # 或者用 -DCMAKE_PREFIX_PATH
cmake --preset msvc2022                               # 或 msvc2019 / msvc2022-v142
cmake --build build\msvc2022 --config Debug
```

### 用 Qt Creator 打开

1. `文件 → 打开文件或项目`，选择根目录的 `CMakeLists.txt`；
2. 选择一个 **Qt 5.15.2 MSVC2019 64bit** 的 Kit（工具链自动选 VS2019/2022）；
3. 构建并运行。若 Kit 配置里已带 Qt 路径，无需任何额外设置（CMakeLists 找不到 Qt 时
   还会自动回退到 `C:/Qt/5.15.2/msvc2019_64`）。

## 目录结构

```
GAMES102/
├── CMakeLists.txt          # CMake 工程（Qt5 + WIN32 + /utf-8）
├── CMakePresets.json       # VS2019 / VS2022（含 v142 工具集）预设
├── configure.bat           # 一键 CMake 配置（VS2022 生成器 + v142 工具集）
├── scripts/
│   ├── build.ps1           # 一键构建（自动探测 VS/Qt + windeployqt）
│   └── run.ps1             # 一键运行
├── src/
│   ├── main.cpp
│   ├── MainWindow.h/.cpp   # 主窗口 + 按钮面板（含 RBF 的 QProcess 调用链）
│   ├── PlotWidget.h/.cpp   # 画布：左键加点、右键删点、绘制曲线
│   └── math/
│       ├── CurveMath.h/.cpp  # 数学库（Lagrange 纯 C++ + 最小二乘 Eigen），与 Qt 解耦
├── tests/
│   └── curve_math_test.cpp   # 数学库自测（无 Qt 依赖，控制台程序）
├── third_party/
│   └── eigen/                # Eigen 3.4.0 头文件（随仓库分发，无需安装）
└── python/                 # uv 环境 + RBF 脚本（pyproject.toml / .venv / src/infer.py）
```

## 扩展新曲线算法

算法都在 `src/math/CurveMath.h`，是**纯 C++、不依赖 Qt** 的函数
（输入 `std::vector<curve::Point>`，输出采样点），方便以后与 Python 结果对拍。
凡是涉及线性代数的部分统一用 **Eigen**（`third_party/eigen`，
`#include <Eigen/Dense>`）。

### 最小二乘类算法（推荐直接用 Eigen）

Eigen 把最容易写错的部分——解方程组——变成一两行。例如已有的最小二乘拟合
（`src/math/CurveMath.cpp`）：

```cpp
Eigen::MatrixXd A(n, m);              // 设计矩阵
Eigen::VectorXd y(n);                 // 目标值
// ... 填充 A(i,j) = x_i^j，y(i) = 数据点纵坐标 ...
Eigen::VectorXd c = A.colPivHouseholderQr().solve(y);   // 解 min ||Ac - y||
```

以后要写的 B 样条拟合、最小二乘曲线、PCA 等算法都是同一套路：

- 解线性方程组：`A.colPivHouseholderQr().solve(b)` / `A.fullPivLu().solve(b)`
- 更稳的最小二乘：`A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b)`
- 教程：https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

### 新增一个算法的标准步骤（已通用化）

现在只需两步：

1. **数学层**：在 `CurveMath.h/.cpp` 加函数（输入 `std::vector<curve::Point>`，
   输出采样点，可自由用 Eigen）；
2. **界面层**：在 `MainWindow` 构造函数里加**一行** `addAlgorithm(...)`。
   幂基插值就是这么加的：

```cpp
addAlgorithm(QStringLiteral("power"), QStringLiteral("插值-幂基"),
             QColor(0x8a, 0x2b, 0xe2),
             [](const std::vector<curve::Point>& pts, int samples) {
                 return curve::powerBaseFunctionInterpolate(pts, samples);
             },
             controls);
```

`addAlgorithm` 自动完成：注册到画布、生成可勾选按钮、「显示所有/隐藏全部」
联动、按各自颜色绘制、自动图例与状态栏提示。**绘图层（PlotWidget）和按钮
逻辑都不需要再改。** 参数化算法（如最小二乘的次数）在 lambda 里通过
`m_plot->approximationDegree()` 读取即可。

3. 可选：在 `tests/curve_math_test.cpp` 补几条断言，跑自测验证数值正确。

### 用自测验证新算法（强烈建议）

`tests/curve_math_test.cpp` 是无 Qt 依赖的控制台自测。加完算法后在这里补几条
断言，跑一遍即可确认数值正确，不用开 GUI 肉眼猜：

```powershell
cmake --build build\msvc2022 --config Debug --target curve_math_test
.\build\msvc2022\Debug\curve_math_test.exe
```

输出 `ALL TESTS PASSED` 即通过。

## Python / RBF 神经网络（作业2）

C++ 通过 QProcess 调用 `uv run --project python python/src/infer.py`：
红点 → `build/io/input.json` → Python 训练 RBF 网络 → `build/io/output.json`
→ 曲线回填到画布（棕色「拟合-RBF神经网络」按钮，隐层数/训练轮数可调）。
**训练过程实时预览**：epoch=0（随机权重）起每 100 轮更新一次曲线，Qt 每
200ms 轮询 `output.json` 刷新，可看到函数从随机逐步拟合到数据点的全过程。
环境配置与手动测试见 [`python/README.md`](python/README.md)。
