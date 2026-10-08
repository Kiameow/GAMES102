# AGENT.md — GAMES102 曲线插值与逼近工具

面向后续开发者 / AI 助手的项目导航。先读本文件，再读 `README.md`（面向用户的安装与使用）。

## 1. 项目是什么

GAMES102 课程作业的**交互式曲线工具**：在画布上用鼠标点出数据点（型值点），一键叠加显示
多种插值 / 逼近 / 参数曲线 / 三次样条算法结果；支持拖动型值点、编辑节点切线（G¹/G⁰）；
另有一个 Python(PyTorch) 实现的 RBF 神经网络拟合（异步回填）。

## 2. 技术栈与构建

| 项 | 说明 |
| --- | --- |
| GUI | Qt 5.15.2（`C:\Qt\5.15.2\msvc2019_64`，MSVC2019 编译） |
| 构建 | CMake ≥ 3.21 + Visual Studio（本机 VS2026 / 18.9，v142 工具集） |
| 线性代数 | Eigen 3.4.0，随仓库分发在 `third_party/eigen/`（纯头文件，无需安装） |
| Python | uv 管理（`python/.venv`），`python/src/infer.py` 做 RBF 训练（PyTorch） |
| 语言 | C++17；源码含中文，MSVC 下统一 `/utf-8` 编译 |

**常用命令**（项目根目录）：

```bat
configure2026.bat                 :: 一键 CMake 配置（VS2026 + v142 工具集，见 CMakePresets.json）
cmake --build --preset msvc2026-v142 --config Debug        :: 构建
set PATH=C:\Qt\5.15.2\msvc2019_64\bin;%PATH%
build\msvc2026-v142\Debug\GAMES102.exe                     :: 运行
cmake --build build\msvc2026-v142 --config Debug --target curve_math_test
build\msvc2026-v142\Debug\curve_math_test.exe              :: 数学库自测（ALL TESTS PASSED）
```

预设（`CMakePresets.json`）：configure 有 `msvc2022 / msvc2022-v142 / msvc2026-v142 / msvc2019`；
build 预设每工具集一条（无固定 configuration，`--config Debug|Release` 构建时指定）。

## 3. 目录结构

```
GAMES102/
├── CMakeLists.txt          # 两个目标：GAMES102(GUI) + curve_math_test(纯C++自测)
├── CMakePresets.json       # VS 各版本/工具集预设
├── configure.bat / configure2026.bat   # 一键配置（VS2022 / VS2026+v142）
├── scripts/build.ps1, run.ps1
├── src/
│   ├── main.cpp
│   ├── MainWindow.h/.cpp   # 窗口 + 控件 + 算法注册（addAlgorithm）
│   ├── PlotWidget.h/.cpp   # 画布：加点/拖点/选点/切线柄、曲线绘制、图例
│   └── math/CurveMath.h/.cpp  # 数学库：所有算法，纯 C++、无 Qt
├── tests/curve_math_test.cpp   # 数学库自测（CHECK 宏，无 Qt 依赖）
├── third_party/eigen/          # Eigen 3.4.0
├── python/                     # uv 工程：src/infer.py（RBF 训练，QProcess 调用）
├── screenshots/                # 截图按钮输出目录（gitignore）
└── hw1.md ~ hw4.md             # 课程作业报告
```

## 4. 架构与核心设计

### 4.1 三层结构

1. **数学层 `curve` 命名空间（`CurveMath.h/.cpp`）**：所有算法。**纯 C++、不依赖 Qt**，
   只依赖 `<Eigen/Dense>`。数据结构只有 `curve::Point{x, y}`。
2. **绘制/交互层 `PlotWidget`**：统一管理"算法层"（`CurveLayer`），负责加点、拖点、选点、
   切线编辑、曲线采样重算与绘制、图例。
3. **窗口层 `MainWindow`**：控制面板（旋钮/下拉框/按钮）+ **算法注册入口**。

### 4.2 算法注册式扩展（最重要的模式）

所有算法统一通过 `MainWindow::addAlgorithm` 注册：

```cpp
using CurveCompute = std::function<std::vector<curve::Point>(
    const std::vector<curve::Point>& pts, int samples)>;

// PlotWidget::registerLayer 登记一层（key/名字/颜色/计算函数），并生成可勾选按钮
addAlgorithm(key, 显示名, QColor, lambda, controlsLayout);
```

- `key`：唯一标识（按钮、状态栏、图层定位用），如 `"lagrange"`、`"spline_natural"`；
- `lambda`：`(pts, samples) -> 采样点向量`——**纯计算，不许碰 UI**；
- 勾选按钮 = 显示该层；`PlotWidget` 在数据/参数变化时自动对**可见**层调用 compute 重算；
- 外部结果层（Python RBF）用 `external=true` 标记，数据变化时由 `clearExternalResults()` 作废。

### 4.3 参数状态放 PlotWidget，算法通过 getter 读取

需要参数的算法，参数状态放在 `PlotWidget`（如 `m_degree / m_sigma / m_lambda /
m_gaussCenters / m_paramBasis`），配套 getter + setter；setter 内部调 `rebuildCurve()` 触发重算。
`MainWindow` 的旋钮信号 → 槽 → `m_plot->setXxx(...)`。算法 lambda 里用 `[this]` 捕获并读
`m_plot->approximationDegree()` 等。

### 4.4 数学库约定

- 函数输入 `std::vector<curve::Point>`（有序），输出**采样点向量**（直接可画）；
- **退化输入一律返回空**：点数 < 2、x 不严格递增（样条）等；
- 解线性方程组统一用 Eigen：最小二乘 `colPivHouseholderQr`，样条三对角 `lu`（钉切线后未必对称）；
- 高次多项式插值（幂基/范德蒙德）会数值病态——新算法别走这条路（用样条/低次/Gauss基）；
- 参数型算法拆 `(t,x)` `(t,y)` 两个标量系统分别拟合再合成，t 来自 `parameterize*`（严格递增）。

### 4.5 交互模型（PlotWidget 鼠标手势）

| 手势 | 行为 |
| --- | --- |
| 左键 | 空处加点；点上有 **500ms 长按** → 拖动型值点（蓝描边） |
| 右键 | 选中最近节点（绿描边）并出现切线柄；空处右键取消选中 |
| 左键拖切线柄端点 | 编辑该节点切线（实时重算曲线） |
| 左键拖靛蓝方块 | 拖动分段 Bezier 层的中间控制点 B1/B2（立即拖，无 500ms；**无右键切线控制**；拖过即 pinned，pts 变化后保留原位） |
| 「删除最近点」按钮 | 删除最近放置的点（右键删点已移除） |

**顶点模式**（MainWindow「顶点模式」下拉框，作用于选中节点）：
- 平滑 = 共享切线（该点 **C¹**，C² 在该点断开）；
- 直线 = 左右共线、长度可独立（**G¹**，共线由 UI `straightenLeft/Right` 保证）；
- 角部 = 左右完全独立（**G⁰**）；
- 自由 = 恢复默认（**C²**，切线由系统解出）。

数学基础：三次样条 = **Hermite 分段插值 + 三转角法**（`hermiteEvaluate` + `cubicSplineNodeSlopes`，
每节点存左右一阶导 `NodeSlope{free, mLeft, mRight}`，自由节点左右相等）。编辑节点 = 钉值，
其 C² 方程被移除、钉值进右端项——**系统自动分块**（一个矩阵一次 `lu().solve()` 求解全部）。

### 4.6 Python RBF 数据协议

`MainWindow` 通过 `QProcess` 启动 `uv run --project python python/src/infer.py`：
红点 → `build/io/input.json`（task/points/params）→ 训练 → `build/io/output.json`
（ok/curve/loss/epoch/done/msg）→ Qt 每 200ms 轮询回填曲线。`build/io/` 是运行期产物。

## 5. 如何添加一个新算法（标准步骤）

1. **数学层**：在 `CurveMath.h/.cpp` 加函数：
   ```cpp
   // 输入有序点列，输出 samples 个采样点；退化输入返回空 vector
   std::vector<Point> myNewAlgorithm(const std::vector<Point>& pts, int samples = 240);
   ```
2. **界面层**：在 `MainWindow` 构造函数加**一行**：
   ```cpp
   addAlgorithm(QStringLiteral("my_key"), QStringLiteral("显示名-说明"),
                QColor(0x00, 0x00, 0x00),  // 颜色避开数据点红色系，见下方色板
                [](const std::vector<curve::Point>& pts, int samples) {
                    return curve::myNewAlgorithm(pts, samples);
                },
                controls);
   ```
3. **要参数**：往 `PlotWidget` 加状态（getter + setter，setter 内 `rebuildCurve()`），
   lambda 改 `[this]` 捕获并从 `m_plot` 读取。
4. **自测**：在 `tests/curve_math_test.cpp` 补断言（CHECK 宏），跑 `curve_math_test.exe`。
5. 构建运行验证。**不需要**改 PlotWidget 绘制逻辑或按钮逻辑——`addAlgorithm` 全包了。

参考现有最小实现：`插值-Lagrange`（无参数）或 `逼近-幂函数最小二乘`（读 degree）。

## 6. 现有算法清单（`MainWindow.cpp` 注册顺序）

| key | 按钮名 | 颜色 | 数学入口 |
| --- | --- | --- | --- |
| `lagrange` | 插值-Lagrange | 蓝 #1f77b4 | `lagrangeInterpolate` |
| `power` | 插值-幂基 | 橙 #ff7f0e | `powerBasePolynomialInterpolate` |
| `gauss` | 插值-Gauss基 | 紫 #9467bd | `gaussBasePolynomialInterpolate`（σ 旋钮） |
| `gauss_ls` | 逼近-Gauss基最小二乘 | 深橄榄 #556b2f | `gaussLeastSquares`（σ + 中心数，k-means 撒点） |
| `least_squares` | 逼近-幂函数最小二乘 | 亮绿 #2ca02c | `leastSquaresPolynomial`（degree） |
| `ridge` | 逼近-岭回归 | 青 #17becf | `ridgeRegression`（degree + λ） |
| `param_uniform` | 参数曲线-均匀参数化 | 粉 #e377c2 | 均匀 t + 基函数选择（幂基/Gauss基） |
| `param_chord` | 参数曲线-弦长参数化 | 灰 #7f7f7f | 弦长 t + 基函数选择 |
| `param_centripetal` | 参数曲线-中心参数化 | 暗金 #b8860b | 中心 t + 基函数选择 |
| `param_foley` | 参数曲线-Foley参数化 | 深紫 #8e44ad | Foley t + 基函数选择 |
| `spline_natural` | 插值-三次样条(自然) | 绿松石 #16a085 | `fitParametricSpline`(Natural) + 切线控制 |
| `spline_clamped` | 插值-三次样条(夹持) | 焦橙 #d35400 | `fitParametricSpline`(Clamped) + 切线控制 |
| `bezier_catmull_rom` | 插值-分段Bezier(Catmull-Rom) | 靛蓝 #4b0082 | `bezierCatmullRomInterpolate`（`buildBezierSegment` + `evaluateBezierSegment`） |
| `chaikin` | 逼近-Chaikin细分(逼近型) | 森林绿 #228b22 | `chaikinSubdivisionCurve`（`chaikinSubdivideOnce` 迭代 4 轮） |
| `bspline_subdiv` | 逼近-均匀三次B样条细分 | 暗青 #008b8b | `bsplineSubdivisionCurve`（`bsplineSubdivideOnce` 迭代 4 轮） |
| `fourpoint_subdiv` | 插值-四点细分 | 灰青 #508b8b | `fourPointSubdivisionCurve`（`fourPointSubdivideOnce` 迭代 4 轮） |
| `rbf`（外部） | 拟合-RBF神经网络(Python) | 棕 #8c564b | Python infer.py 异步回填 |

参数型四种按钮共享「参数拟合基」下拉框（幂基/ Gauss基），四次样条按钮共享「顶点模式」。

## 7. 关键约定与注意事项

- **曲线颜色**：数据点是红色（#d62728），曲线慎用红色系；新增算法参考色板
  蓝/橙/绿/紫/青/棕/粉/灰（已用）+ 深橄榄/暗金/深紫/绿松石/焦橙（补充色）。
- **configure.bat 必须纯 ASCII + CRLF**（cmd 解析要求），编辑时保持。
- **中文源码**：MSVC 用 `/utf-8`（CMakeLists 已加），新文件保持 UTF-8。
- **截图**：控制面板「保存截图(PNG)」→ `screenshots/GAMES102_yyyy-MM-dd_HHmmss.png`
  （时间戳命名，路径复制到剪贴板）；`screenshots/` 已 gitignore。
- **WIP / 已知问题**：
  - 分段 Bezier（`bezierCatmullRomInterpolate`）已实现并接入 GUI：每段中间控制点
    `BezierSegment{Point control[2]}` 由 Catmull-Rom 切线构造，边界段用镜像虚拟点
    （`catmullRomVirtualPoint`），n==2 时直接退化为直线段；早期的
    `evaluateCubicBernsteinBasisCurve` 笔误（`pow(1-t, 4-i)` → `3-i`）已修复。
    控制点由 PlotWidget 全局管理（`m_bezierSegments` + `m_bezierPinned`，与 pts 平行）：
    勾选该层时显示靛蓝方块 + 虚线控制多边形，左键可直接拖动（拖过即 pinned，
    pts 变化后保留原位），曲线用 `bezierCatmullRomInterpolate(pts, segs, samples)`
    读取拖过的控制点。
  - Chaikin 细分（`chaikinSubdivisionCurve`，逼近型）已接入 GUI：每轮每条边生成两个
    割角点（3/4:1/4、1/4:3/4），老点全部抛弃，开放曲线保留首尾端点；默认迭代 4 轮，
    直接把最终多边形画成逼近曲线，不显示中间细分点（极限为二次 B 样条）。
  - 均匀三次 B 样条细分（`bsplineSubdivisionCurve`，逼近型）已接入 GUI：每轮内部顶点
    1/8:3/4:1/8 平滑 + 每条边 1/2:1/2 中点（ν′_{2i}=1/8ν_{i−1}+3/4ν_i+1/8ν_{i+1}、
    ν′_{2i+1}=1/2ν_i+1/2ν_{i+1}），老点抛弃、保留首尾端点（钳制端），极限为均匀三次
    B 样条。与 Chaikin 共用细分迭代骨架 `subdivideIterate`（函数指针传入单轮规则）。
  - 4 点插值细分（`fourPointSubdivisionCurve`，插值型，Dyn–Levin–Gregory w=1/16）已接入
    GUI：每轮旧点全部保留、每条边中点插入一个新点（9/16:9/16:−1/16:−1/16 模板 = 过 4
    个相邻点的三次多项式在中点的值），边界用镜像虚拟点补齐（n==3 时两条边都是边界边，
    各用一个虚拟点）；极限 C¹ 且过全部原始点，能精确还原三次多项式。迭代骨架同样复用
    `subdivideIterate`。
  - **增量构建坑**：改了 `.h` 后 MSBuild 偶尔不重编依赖它的 .cpp（obj 比源码新）。若链接报
    `unresolved external symbol`，先查 `build\msvc2026-v142\<Target>.dir\Debug\*.obj` 时间戳，
    删掉对应 obj 或 touch 源码强制重编译。
  - 高次幂基插值（拉满次数）数值病态会"乱飞"，属预期，不是 bug。
- **作业报告**：`hw1.md ~ hw4.md`，结论要点见各文件。

## 8. 测试

`tests/curve_math_test.cpp` 用 `CHECK(cond)` 宏（失败计数），目标 `curve_math_test` 无 Qt
依赖，跑 `ALL TESTS PASSED` 即过。**新增算法强烈建议补断言**（精确恢复已知函数 / 退化输入
返回空 / 无 NaN）。
