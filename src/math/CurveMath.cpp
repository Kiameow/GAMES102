#include "CurveMath.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace curve {

namespace {
constexpr double kEps = 1e-12;
}

// ---------------------------------------------------------------------------
// Lagrange 插值（纯 C++，无需 Eigen）
// ---------------------------------------------------------------------------

double lagrangeAt(const std::vector<Point>& pts, double x) {
    const int n = static_cast<int>(pts.size());
    if (n == 0) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        // 采样点恰好落在数据点上时直接返回，避免 0/0
        if (std::abs(x - pts[i].x) < kEps) return pts[i].y;
        double term = pts[i].y;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            term *= (x - pts[j].x) / (pts[i].x - pts[j].x);
        }
        sum += term;
    }
    return sum;
}

std::vector<Point> lagrangeInterpolate(const std::vector<Point>& pts, int samples) {
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    double xmin = pts.front().x;
    double xmax = xmin;
    for (const Point& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
    }
    if (xmax - xmin < kEps) return out;  // 所有 x 相同，无法插值

    samples = std::max(samples, 2);
    out.reserve(static_cast<std::size_t>(samples));
    for (int k = 0; k < samples; ++k) {
        const double x = xmin + (xmax - xmin) * static_cast<double>(k) / (samples - 1);
        out.push_back({x, lagrangeAt(pts, x)});
    }
    return out;
}

// 基函数插值法和逼近法都需要求解待定系数，而拉格朗日插值法则不需要
// ---------------------------------------------------------------------------
// 幂基函数插值，可以复用最小二乘法
// ---------------------------------------------------------------------------
std::vector<Point> powerBasePolynomialInterpolate(const std::vector<Point>& pts, int samples) {
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    int degree = pts.size() - 1;
    const std::vector<double> coefs = fitPolynomial(pts, degree);
    if (coefs.empty()) return out;

    double xmin = pts.front().x;
    double xmax = xmin;
    for (const Point& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
    }

    samples = std::max(samples, 2);
    out.reserve(static_cast<std::size_t>(samples));
    for (int k = 0; k < samples; ++k) {
        const double x = xmin + (xmax - xmin) * static_cast<double>(k) / (samples - 1);
        out.push_back({x, evaluatePolynomial(coefs, x)});
    }
    return out;
}

// ---------------------------------------------------------------------------
// Gauss基函数插值
// ---------------------------------------------------------------------------
double gaussNodeValue(double insertPointX, double nodeX, double sigma) {
    return std::exp(-(insertPointX - nodeX) * (insertPointX - nodeX) / (2 * sigma * sigma));
}

// 求解高斯基函数插值的待定系数。
// 未知数：n 个高斯系数 cⱼ（每个数据点一个中心）+ 1 个常数 c₀ = n+1 个；
// 方程：n 个插值条件 f(xᵢ) = yᵢ + 1 个约束 Σcⱼ = 0。
// 增广系统：
//     [ Φ   1 ] [c ]   [y]
//     [ 1ᵀ  0 ] [c₀] = [0]
// 返回系数 [c₀, ..., c_{n-1}, 常数项]，长度 n+1。
std::vector<double> fitGaussBasePolynomial(const std::vector<Point>& pts, double sigma) {
    std::vector<double> coefs;
    const int n = static_cast<int>(pts.size());
    if (n < 2 || sigma <= 0.0) return coefs;

    const int m = n + 1;  // 未知数个数
    Eigen::MatrixXd A(m, m);
    Eigen::VectorXd b(m);

    // 前 n 行：插值条件。Φ(i,j) = φ(xᵢ − xⱼ)，最后一列是常数项
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) A(i, j) = gaussNodeValue(pts[i].x, pts[j].x, sigma);
        A(i, n) = 1.0;
        b(i) = pts[i].y;
    }
    // 最后一行：约束 Σcⱼ = 0，保证增广系统非奇异、解唯一
    for (int j = 0; j < n; ++j) A(n, j) = 1.0;
    A(n, n) = 0.0;
    b(n) = 0.0;

    const Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

    coefs.reserve(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) coefs.push_back(c(j));
    return coefs;
}

// coefs[0..n-1] = 各高斯基的系数，coefs[n] = 常数项
double evaluateGaussBasePolynomial(const std::vector<double>& coefs,
                                   const std::vector<Point>& pts, double sigma, double x) {
    const int n = static_cast<int>(pts.size());
    if (coefs.size() < static_cast<std::size_t>(n) + 1) return 0.0;
    double y = coefs[n];
    for (int i = 0; i < n; ++i) y += coefs[i] * gaussNodeValue(x, pts[i].x, sigma);
    return y;
}

std::vector<Point> gaussBasePolynomialInterpolate(const std::vector<Point>& pts, double sigma, int samples) {
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    const std::vector<double> coefs = fitGaussBasePolynomial(pts, sigma);
    if (coefs.empty()) return out;

    double xmin = pts.front().x;
    double xmax = xmin;
    for (const Point& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
    }

    samples = std::max(samples, 2);
    out.reserve(static_cast<std::size_t>(samples));
    for (int k = 0; k < samples; ++k) {
        const double x = xmin + (xmax - xmin) * static_cast<double>(k) / (samples - 1);
        out.push_back({ x, evaluateGaussBasePolynomial(coefs, pts, sigma, x) });
    }
    return out;
}

// ---------------------------------------------------------------------------
// 最小二乘多项式逼近（Eigen 实现）
// ---------------------------------------------------------------------------

std::vector<double> fitPolynomial(const std::vector<Point>& pts, int degree) {
    std::vector<double> coefs;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return coefs;

    degree = std::clamp(degree, 1, n - 1);  // 次数限制在 [1, n-1]
    const int m = degree + 1;               // 未知数个数

    // 设计矩阵 A(i,j) = x_i^j，右端项 y_i
    Eigen::MatrixXd A(n, m);
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        A(i, 0) = 1.0;
        for (int j = 1; j < m; ++j) A(i, j) = A(i, j - 1) * pts[i].x; // 前项累乘减少计算量
        y(i) = pts[i].y;
    }

    // 求解最小二乘问题 min || A c - y ||：列主元 QR 分解（数值稳定）
    // 最小二乘法得到的目标函数可以简化为 || Ax ||2 - 2 * transpose(Ax) * b + transpose(b) * b
    // 求导之后，令等于0，得到 transpose(A) * A * x = transpose(A) * b
    // 如果对A作QR经济型分解（要求A列满秩），则求导之后令等于0可以得到transpose(R) * R * x = transpose(R) * transpose(Q) * b
    // QR经济型分解可以得到一个可逆的R，以及Q * transpose(Q) = I(m*m，但是秩只有n)，因此可以变为Rx = transpose(Q)*b
    // 同时左乘Q之后，得到Ax = Pb，其中P = Q * transpose(Q)

    // 这里colPivHouseholderQr会智能检测A的长宽，如果是瘦长型的，那就会做经济型分解，即上述的最小而成法过程，而不是插值过程
    const Eigen::VectorXd c = A.colPivHouseholderQr().solve(y);

    coefs.reserve(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) coefs.push_back(c(j));
    return coefs;
}

double evaluatePolynomial(const std::vector<double>& coefs, double x) {
    double y = 0.0;
    // Horner 法则：从最高次项往低次累乘
    for (auto it = coefs.rbegin(); it != coefs.rend(); ++it) y = y * x + *it;
    return y;
}

std::vector<Point> leastSquaresPolynomial(const std::vector<Point>& pts,
                                          int degree,
                                          int samples) {
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    const std::vector<double> coefs = fitPolynomial(pts, degree);
    if (coefs.empty()) return out;

    double xmin = pts.front().x;
    double xmax = xmin;
    for (const Point& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
    }

    samples = std::max(samples, 2);
    out.reserve(static_cast<std::size_t>(samples));
    for (int k = 0; k < samples; ++k) {
        const double x = xmin + (xmax - xmin) * static_cast<double>(k) / (samples - 1);
        out.push_back({x, evaluatePolynomial(coefs, x)});
    }
    return out;
}

// ---------------------------------------------------------------------------
// 岭回归（L2 正则最小二乘）：min ||Ac - y||² + λ·||c||²，把第二项尝试合进第一项，把A变成一个增广矩阵，
// 用QR解就可以
// ---------------------------------------------------------------------------

std::vector<double> fitPolynomialRidge(const std::vector<Point>& pts, int degree, double lambda) {
    std::vector<double> coefs;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return coefs;

    degree = std::clamp(degree, 1, n - 1);  // 次数限制在 [1, n-1]
    const int m = degree + 1;               // 未知数个数
    const double s = std::sqrt(std::max(0.0, lambda));  // 增广块里的 √λ

    // 增广矩阵：n 个数据行 + (m-1) 个正则行（跳过常数项 c₀，不惩罚截距）
    const int rows = n + (m - 1);
    Eigen::MatrixXd A(rows, m);
    Eigen::VectorXd b(rows);

    // 数据行（与 fitPolynomial 相同）
    for (int i = 0; i < n; ++i) {
        A(i, 0) = 1.0;
        for (int j = 1; j < m; ++j) A(i, j) = A(i, j - 1) * pts[i].x; // 前项累乘减少计算量
        b(i) = pts[i].y;
    }
    // 正则行：第 k 行只在第 (k+1) 列放 √λ，其余为 0，右端补 0
    // 若想连常数项一起惩罚：把 (m-1) 改成 m，对角线放到第 k 列
    for (int k = 0; k < m - 1; ++k) {
        A(n + k, 0) = 0.0;
        for (int j = 1; j < m; ++j) A(n + k, j) = (j == k + 1) ? s : 0.0;
        b(n + k) = 0.0;
    }

    const Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

    coefs.reserve(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) coefs.push_back(c(j));
    return coefs;
}

double evaluatePolynomialRidge(const std::vector<double>& coefs, double x) {
    double y = 0.0;
    // Horner 法则：从最高次项往低次累乘
    for (auto it = coefs.rbegin(); it != coefs.rend(); ++it) y = y * x + *it;
    return y;
}

std::vector<Point> ridgeRegression(const std::vector<Point>& pts,
                                   int degree,
                                   double lambda,
                                   int samples) {
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    const std::vector<double> coefs = fitPolynomialRidge(pts, degree, lambda);
    if (coefs.empty()) return out;

    double xmin = pts.front().x;
    double xmax = xmin;
    for (const Point& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
    }

    samples = std::max(samples, 2);
    out.reserve(static_cast<std::size_t>(samples));
    for (int k = 0; k < samples; ++k) {
        const double x = xmin + (xmax - xmin) * static_cast<double>(k) / (samples - 1);
        out.push_back({x, evaluatePolynomialRidge(coefs, x)});
    }
    return out;
}

}  // namespace curve
