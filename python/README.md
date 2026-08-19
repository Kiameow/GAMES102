# python/ — RBF 神经网络拟合（作业2）

Qt(C++) 负责采集数据点并绘制，Python 负责训练 RBF 神经网络，二者通过
**QProcess + JSON** 通信（协议见 `python/src/infer.py` 顶部注释）。

## 环境

```powershell
cd python
uv sync          # 用 uv 按 pyproject.toml 创建 .venv 并安装 numpy/torch(cpu)
uv run python -c "import torch; print(torch.__version__)"   # 验证
```

## 手动测试（不经过 Qt）

```powershell
# 在项目根目录执行（脚本在 python/src/ 下，环境在 python/ 下）
uv run --project python python/src/infer.py build/io/input.json build/io/output.json
```

`infer.py` 逻辑：读 input.json → 构建 RBF 网络（centers/σ/权重全部可学习，
Adam + MSE 端到端训练）→ 在数据范围内采样 240 点 → 写 output.json。

**实时预览**：训练一开始（epoch=0，随机权重）就写一版 output.json，
之后每 `report_interval`（默认 100）轮更新一次（`done:false`），训练结束写
`done:true` 的最终结果。Qt 侧每 200ms 轮询一次，所以你能看到曲线从随机
函数逐步拟合到数据点的全过程。

## 常见问题

- **torch 装不上 / 下载太大**：pyproject.toml 里配了清华镜像 + CPU wheel 源，
  若失效换官方源 `https://download.pytorch.org/whl/cpu`
- **训练发散（NaN）**：调小 lr（0.001）或调大 sigma_init（0.5）
- **曲线太平**：隐层数太少或 σ 太大；**曲线过冲**：σ 太小或 epochs 太多
