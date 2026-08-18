#pragma once

#include <vector>

namespace curve {

// 二维数据点（纯数值结构，与 Qt 无关）
struct Point {
    double x = 0.0;
    double y = 0.0;
};

// ---------------------------------------------------------------------------
// 多项式插值（Lagrange 形式）
// 纯 C++ 实现：插值本质是"基函数加权求和"，不需要线性代数，
// 因此不引入 Eigen（Eigen 用在需要解方程组/最小二乘的地方，见下方）。
// ---------------------------------------------------------------------------

// 在 x 处求 Lagrange 插值多项式的值
double lagrangeAt(const std::vector<Point>& pts, double x);

// Lagrange 插值：返回经过所有 pts 的曲线，共 samples 个采样点。
// 少于 2 个点或所有 x 相同时返回空。
std::vector<Point> lagrangeInterpolate(const std::vector<Point>& pts, int samples = 240);

std::vector<Point> powerBasePolynomialInterpolate(const std::vector<Point>& pts, int samples = 240);

// Gauss 基函数插值（RBF 插值）：
//     f(x) = Σⱼ cⱼ·φ(x − xⱼ) + c₀，φ(r) = exp(−r²/(2σ²))
// 每个数据点一个高斯中心（n 个基），加一个常数项 c₀，
// 增广系统 [Φ 1; 1ᵀ 0] 求解（约束 Σcⱼ = 0 保证解唯一），
// 曲线精确经过全部 n 个点。σ 是高斯宽度，控制曲线的"局部性"。
double gaussNodeValue(double insertPointX, double nodeX, double sigma);
std::vector<double> fitGaussBasePolynomial(const std::vector<Point>& pts, double sigma);
double evaluateGaussBasePolynomial(const std::vector<double>& coefs, const std::vector<Point>& pts, double sigma, double x);
std::vector<Point> gaussBasePolynomialInterpolate(const std::vector<Point>& pts, double sigma, int samples = 240);

// ---------------------------------------------------------------------------
// 最小二乘多项式逼近（基于 Eigen）
// 底层用 Eigen::MatrixXd 组装设计矩阵，Eigen::ColPivHouseholderQR 求解，
// 数值稳定：重复 x、病态点都不会崩溃（返回最小二乘意义下的解）。
// ---------------------------------------------------------------------------

// 最小二乘多项式拟合：用 degree 次多项式
//     y = c[0] + c[1]*x + c[2]*x^2 + ... + c[degree]*x^degree
// 拟合 pts，返回系数向量 c（长度 = degree+1）。
// degree 会自动限制在 [1, n-1]；少于 2 个点时返回空。
std::vector<double> fitPolynomial(const std::vector<Point>& pts, int degree);

double evaluatePolynomial(const std::vector<double>& coefs, double x);

// 最小二乘多项式逼近曲线：对 fitPolynomial 得到的多项式采样，
// 返回 samples 个点（直接给 PlotWidget 绘制用）。
std::vector<Point> leastSquaresPolynomial(const std::vector<Point>& pts,
                                          int degree,
                                          int samples = 240);

// ---------------------------------------------------------------------------
// 岭回归（L2 正则最小二乘）：min ||Ac - y||² + λ·||c||²
// 实现：增广矩阵法 [A; √λ·D] c ≈ [y; 0]，直接走 QR 求解，数值稳定
//       （比法方程 (AᵀA + λI)c = Aᵀy 不平方条件数）。
// 默认不惩罚常数项 c₀（截距），只惩罚 c₁..c_{m-1}；λ=0 退化为最小二乘。
// ---------------------------------------------------------------------------
std::vector<double> fitPolynomialRidge(const std::vector<Point>& pts, int degree, double lambda);

double evaluatePolynomialRidge(const std::vector<double>& coefs, double x);

std::vector<Point> ridgeRegression(const std::vector<Point>& pts,
                                   int degree,
                                   double lambda,
                                   int samples = 240);

}  // namespace curve
