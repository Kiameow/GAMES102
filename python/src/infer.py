#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RBF 神经网络拟合 —— 供 C++(Qt) 通过 QProcess 调用。

数据协议（JSON，全部 UTF-8）：
  input.json:
    {
      "task": "fit_rbf",
      "points": [[x0,y0],[x1,y1],...],
      "params": {"n_centers":6, "sigma_init":0.3, "epochs":2000, "lr":0.01,
                 "samples":240, "seed":0, "report_interval":100}
    }
  output.json（训练中会不断更新；写入是原子的，Qt 轮询读不会读到半个文件）:
    {"ok": true, "curve": [[x,y],...], "loss": 0.0004, "epoch": 100,
     "done": false, "msg": "epoch=100 loss=0.000412"}     <- 中间进度
    {"ok": true, "curve": [[x,y],...], "loss": 0.0001, "epoch": 2000,
     "done": true,  "msg": "..."}                        <- 训练完成
    {"ok": false, "msg": "错误信息"}                      <- 失败（exit code = 1）

进度报告：epoch=0（随机初始化）先出一版，之后每 report_interval 轮更新一次，
让 Qt 侧能实时看到拟合过程。用法（在项目根目录执行）：
  uv run --project python python/src/infer.py [input.json] [output.json]
"""
import json
import os
import sys

import numpy as np
import torch
import torch.nn as nn


def find_project_root() -> str:
    """向上查找项目根（以 CMakeLists.txt 为标记），脚本位置变动也不受影响。"""
    d = os.path.dirname(os.path.abspath(__file__))
    while True:
        if os.path.exists(os.path.join(d, "CMakeLists.txt")):
            return d
        parent = os.path.dirname(d)
        if parent == d:  # 已到盘符根仍未找到，退回当前目录
            return d
        d = parent


PROJECT_ROOT = find_project_root()


def resolve(path: str) -> str:
    """相对路径按项目根目录解析，绝对路径原样返回。"""
    return os.path.join(PROJECT_ROOT, path) if not os.path.isabs(path) else path


def write_output(path: str, ok: bool, **fields) -> None:
    """写 output.json。
    注意：不用"临时文件 + os.replace"。Windows 下 os.replace 会与 Qt 侧轮询
    打开的文件（未带 FILE_SHARE_DELETE）冲突抛 PermissionError；而 Qt 轮询
    读到写了一半的 JSON 会解析失败自动跳过，所以直接截断写入即可。"""
    payload = {"ok": ok, **fields}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


class RBFNet(nn.Module):
    """RBF 网络：f(x) = Σⱼ wⱼ·φ(||x − cⱼ||)，φ(r) = exp(−r²/(2σ²))
    centers / σ / 输出权重 w 全部可学习（端到端梯度下降）。"""

    def __init__(self, centers: np.ndarray, sigma: float):
        super().__init__()
        # 定义在init里的变量都是可以学习的变量，但是如何学习是我后续需要重点理解的内容 ！！
        self.centers = nn.Parameter(torch.from_numpy(centers))       # (n_c, 1)
        self.log_sigma = nn.Parameter(torch.log(torch.tensor(sigma, dtype=torch.float32)))
        # 输出层的定义其实就是声明了一个矩阵运算，像这里就是声明最后一个隐藏层到输出，应该经过一个右乘centers * 1的权重向量
        self.out = nn.Linear(centers.shape[0], 1, bias=True)  

    def forward(self, x: torch.Tensor) -> torch.Tensor:              # x: (N,1)
        # d2的计算用到了广播运算，非常值得学习，[pts, 1, 1] - [1, centers, 1]
        # 最后会变成一个[pts, centers, 1]的矩阵，非常快速的完成了每个数据点和中心的x之差运算
        # .sum(-1)不会改变数据，只是少了最后一维
        d2 = (x[:, None, :] - self.centers[None, :, :]).pow(2).sum(-1)
        sigma = self.log_sigma.exp()
        phi = torch.exp(-d2 / (2.0 * sigma * sigma))
        return self.out(phi)                                         # (N,1)


def fit_rbf(points: np.ndarray, params: dict, report=None):
    """对 points (N,2) 训练 RBF 网络。
    report(epoch, loss, curve) 是进度回调：epoch=0 时用随机初始化权重
    先报告一版，之后每 report_interval 轮报告一次。返回 (最终曲线, 最终损失)。"""
    x = points[:, 0:1].astype(np.float32)
    y = points[:, 1:2].astype(np.float32)

    n_centers = max(1, min(int(params.get("n_centers", 6)), 200))
    sigma_init = max(1e-3, float(params.get("sigma_init", 0.3)))
    epochs = max(1, min(int(params.get("epochs", 2000)), 200000))
    lr = float(params.get("lr", 0.01))
    samples = max(2, min(int(params.get("samples", 240)), 5000))
    report_interval = max(1, int(params.get("report_interval", 100)))
    seed = int(params.get("seed", 0))
    torch.manual_seed(seed)
    np.random.seed(seed)

    # 隐层中心初始化：在 x 范围内等距撒点（点数少时够用；点密集可换 k-means）
    centers = np.linspace(float(x.min()), float(x.max()), n_centers, dtype=np.float32)
    centers = centers.reshape(-1, 1)

    model = RBFNet(centers, sigma_init)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    xt = torch.from_numpy(x)
    yt = torch.from_numpy(y)

    # 采样网格：训练过程中反复用它出中间曲线
    xs = np.linspace(float(x.min()), float(x.max()), samples, dtype=np.float32)

    def sample_curve():
        with torch.no_grad(): # 禁用梯度下降，防止权重在推理时被更新
            ys = model(torch.from_numpy(xs.reshape(-1, 1))).numpy().ravel()
        if not np.all(np.isfinite(ys)):
            raise RuntimeError("训练发散（输出出现 NaN/Inf），请减小 lr 或增大 sigma_init")
        return [[float(a), float(b)] for a, b in zip(xs, ys)]

    # 第 0 轮：随机初始化权重就出一版，让 Qt 立刻能看到函数
    if report:
        report(0, None, sample_curve())

    loss = None
    for epoch in range(1, epochs + 1):
        opt.zero_grad()
        # loss的类型： tensor(0.0042, grad_fn=<MseLossBackward0>)
        # 只有知道了梯度函数的具体构造，才能成功反向传播，知道哪个变量要被更新
        # 输出报告时，可以用loss.detach()来丢弃grad_fn信息
        loss = nn.functional.mse_loss(model(xt), yt) # 此处调用了forward方法，计算目标函数值
        loss.backward()
        opt.step()
        if report and epoch % report_interval == 0:
            report(epoch, float(loss.detach()), sample_curve())

    return sample_curve(), float(loss.detach())


def main() -> None:
    input_path = resolve(sys.argv[1]) if len(sys.argv) > 1 else resolve("build/io/input.json")
    output_path = resolve(sys.argv[2]) if len(sys.argv) > 2 else resolve("build/io/output.json")

    try:
        with open(input_path, "r", encoding="utf-8") as f:
            req = json.load(f)

        task = req.get("task", "fit_rbf")
        if task != "fit_rbf":
            raise RuntimeError(f"未知任务: {task}")

        points = np.asarray(req["points"], dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 2:
            raise RuntimeError("points 必须是 [[x,y],...] 形状")
        if points.shape[0] < 2:
            raise RuntimeError("至少需要 2 个数据点")

        params = req.get("params", {})

        # 训练中的进度回调：写 done=false 的中间结果
        def progress(epoch, loss, curve):
            msg = (f"epoch={epoch} loss={loss:.6g}" if loss is not None
                   else "epoch=0（随机初始化）")
            write_output(output_path, True, curve=curve, loss=loss,
                         epoch=epoch, done=False, msg=msg)

        curve, loss = fit_rbf(points, params, report=progress)

        # 训练完成：写 done=true 的最终结果
        msg = (f"RBF-NN n_centers={params.get('n_centers', 6)} "
               f"epochs={params.get('epochs', 2000)} loss={loss:.6g}")
        write_output(output_path, True, curve=curve, loss=loss,
                     epoch=int(params.get("epochs", 2000)), done=True, msg=msg)
        print(msg)
    except Exception as e:  # 任何错误都写成 output.json 并返回非零退出码
        write_output(output_path, False, msg=f"{type(e).__name__}: {e}")
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
