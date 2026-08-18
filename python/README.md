# Python 神经网络 I/O 规划（暂未实现）

> 本文档只做**方案规划**，方便后续 GAMES102 作业（例如用神经网络做曲线拟合 /
> 点云生成 / 参数化）时参考。当前 Qt 程序不依赖 Python，无需安装任何东西。

## 目标场景

- 用户在 Qt 画布上点出数据点（红点）→ 调用 Python 里的神经网络
  （训练好的模型 / 实时训练）→ 把网络输出的曲线/结果画回 Qt 窗口。
- 典型用例：给定离散点，用训练好的网络输出更平滑的曲线、控制点、
  或预测形状（对应 GAMES102 的曲线/曲面拟合、生成类作业）。

## 推荐方案：QProcess + JSON（简单、可靠、易调试）

```
Qt (C++ 主程序)                    Python (独立进程)
─────────────────                  ─────────────────
1. 用户点好红点
2. 写出 build/io/input.json
        ────────────────────────►
   QProcess 启动 python.exe
   python/infer.py input.json     3. 读取 input.json
                                  4. 加载模型 / 运行网络
                                  5. 写出 build/io/output.json
        ◄────────────────────────
6. 等待进程结束，读取 output.json
7. 把结果曲线画到画布上
```

优点：两个进程完全解耦，协议（JSON）清晰，任何一步出错都能单独跑
`python infer.py` 调试；不引入 pybind11 的编译复杂度。

### 数据协议（草案）

`build/io/input.json`:

```json
{
  "task": "approximate",            // 或 "interpolate" / "generate"
  "points": [[0.1, 0.2], [0.3, 0.8], [0.6, 0.4], [0.9, 0.7]],
  "params": { "degree": 3 }         // 可选参数
}
```

`build/io/output.json`:

```json
{
  "ok": true,
  "curve": [[0.1, 0.21], [0.12, 0.24], "...", [0.9, 0.71]]
}
```

Qt 端用 `QProcess` 调用，示例结构（后续实现）：

```cpp
QProcess py;
py.start("python", {"python/infer.py", "build/io/input.json"});
// finished 信号里读取 build/io/output.json，转成 QPolygonF 画出来
```

### 目录规划

```
python/
├── README.md            # 本文档
├── requirements.txt     # numpy / torch 等（按需）
├── infer.py             # 推理入口：读 input.json -> 跑网络 -> 写 output.json
├── train.py             # 训练脚本（离线跑，不占用 Qt）
├── models/              # 模型定义与权重
└── .venv/               # 虚拟环境（已加入 .gitignore）
```

## 备选方案（按需升级）

| 方案 | 优点 | 缺点 | 适用 |
| --- | --- | --- | --- |
| QProcess + JSON（推荐） | 简单、易调试、协议清晰 | 每帧都有进程开销 | 离线/准实时推理 |
| QProcess 长驻 + 管道 | 免去反复启动 | 协议解析稍复杂 | 需要实时交互 |
| pybind11 嵌入 | 零进程开销、可传张量 | 编译配置复杂、依赖捆绑 | 性能敏感、大量交互 |
| 本地 HTTP (FastAPI) | 跨语言最通用、可远程 | 多一层网络栈 | 多人协作/服务化 |

> 入门阶段用 QProcess + JSON 就足够，等作业明确需要实时训练可视化时再升级。

## 后续实现清单

- [ ] 在 Qt 中把当前红点导出为 `input.json`
- [ ] 写 `python/infer.py` 骨架（读 JSON、跑网络、写 JSON）
- [ ] 在 `MainWindow` 加「调用 Python」按钮 + 结果回显
- [ ] 训练一个最简单的曲线拟合网络（MLP 拟合点云）
- [ ] （可选）进程长驻 + 增量输入
