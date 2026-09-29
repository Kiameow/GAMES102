#include "CurveMath.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace curve {

namespace {
constexpr double kEps = 1e-12;

// 单轮细分函数指针（Chaikin / B 样条共用同一细分迭代骨架）
using SubdivideStepFn = std::vector<Point> (*)(const std::vector<Point>&);

// 通用细分迭代：对 pts 迭代 iterations 轮，每轮用 step 生成新点列，
// 返回最终多边形（直接作为逼近曲线）。退化输入（<2 点）返回空。
std::vector<Point> subdivideIterate(const std::vector<Point>& pts, int iterations,
                                    SubdivideStepFn step) {
    if (pts.size() < 2) return {};
    std::vector<Point> cur = pts;
    iterations = std::clamp(iterations, 0, 10);
    for (int k = 0; k < iterations; ++k) {
        cur = step(cur);
        if (cur.size() < 3) break;  // 只剩直线段，再细分无意义（n==2 时每轮原样返回）
    }
    return cur;
}
}

double L2(const Point& p1, const Point& p2)
{
    const double dx = p2.x - p1.x;
    const double dy = p2.y - p1.y;
    return std::sqrt(dx * dx + dy * dy);
}

double angleInRadians(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2)
{
    double norm1 = v1.norm();
    double norm2 = v2.norm();
    // 退化边（长度≈0，如 Foley 边界处 p0==p1 或 p2==p3）：没有转角信息，
    // 视为"无转角"，返回 0，避免 0/0=NaN 污染后续计算。
    if (norm1 < kEps || norm2 < kEps) return 0.0;

    double dot = v1.dot(v2);
    double cos_angle = dot / (norm1 * norm2);

    cos_angle = std::clamp(cos_angle, -1.0, 1.0);

    return std::acos(cos_angle);
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
// Gauss 基函数最小二乘拟合（逼近，非插值；k-means 撒点）
// ---------------------------------------------------------------------------

std::vector<double> kMeansCenters1D(const std::vector<double>& xs, int k) {
    std::vector<double> centers;
    const int n = static_cast<int>(xs.size());
    if (n == 0 || k <= 0) return centers;
    k = std::min(k, n);

    // 初始化：在 [xmin, xmax] 均匀撒 k 个中心（k=1 时取中点）
    const auto [xminIt, xmaxIt] = std::minmax_element(xs.begin(), xs.end());
    const double xmin = *xminIt;
    const double xmax = *xmaxIt;
    centers.resize(static_cast<std::size_t>(k));
    for (int j = 0; j < k; ++j) {
        centers[j] = xmin + (xmax - xmin) * (k == 1 ? 0.5 : static_cast<double>(j) / (k - 1));
    }

    std::vector<int> assign(static_cast<std::size_t>(n));
    std::vector<double> sum(static_cast<std::size_t>(k));
    std::vector<int> cnt(static_cast<std::size_t>(k));
    for (int iter = 0; iter < 100; ++iter) {
        // E 步：每个点归到最近中心
        for (int i = 0; i < n; ++i) {
            double best = std::numeric_limits<double>::max();
            int bj = 0;
            for (int j = 0; j < k; ++j) {
                const double d = std::abs(xs[static_cast<std::size_t>(i)] - centers[static_cast<std::size_t>(j)]);
                if (d < best) { best = d; bj = j; }
            }
            assign[static_cast<std::size_t>(i)] = bj;
        }
        // M 步：每个簇取均值；空簇保持原中心不动
        std::fill(sum.begin(), sum.end(), 0.0);
        std::fill(cnt.begin(), cnt.end(), 0);
        for (int i = 0; i < n; ++i) {
            sum[static_cast<std::size_t>(assign[static_cast<std::size_t>(i)])] += xs[static_cast<std::size_t>(i)];
            cnt[static_cast<std::size_t>(assign[static_cast<std::size_t>(i)])]++;
        }
        bool changed = false;
        for (int j = 0; j < k; ++j) {
            if (cnt[static_cast<std::size_t>(j)] == 0) continue;
            const double nc = sum[static_cast<std::size_t>(j)] / cnt[static_cast<std::size_t>(j)];
            if (std::abs(nc - centers[static_cast<std::size_t>(j)]) > 1e-12) changed = true;
            centers[static_cast<std::size_t>(j)] = nc;
        }
        if (!changed) break;
    }
    std::sort(centers.begin(), centers.end());
    return centers;
}

std::vector<double> fitGaussLeastSquares(const std::vector<Point>& pts, double sigma,
                                         int numCenters, std::vector<double>& centers) {
    centers.clear();
    std::vector<double> coefs;
    const int n = static_cast<int>(pts.size());
    if (n < 2 || sigma <= 0.0 || numCenters <= 0) return coefs;

    // 1) 用 1D k-means 在 x 坐标上撒点得到中心（自适应数据分布）
    std::vector<double> xs(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) xs[static_cast<std::size_t>(i)] = pts[static_cast<std::size_t>(i)].x;
    centers = kMeansCenters1D(xs, std::min(numCenters, n));
    const int m = static_cast<int>(centers.size());
    if (m < 1) return coefs;

    // 2) 设计矩阵 A(i,j) = φ(xᵢ − centerⱼ)，最小二乘解 c
    Eigen::MatrixXd A(n, m);
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j)
            A(i, j) = gaussNodeValue(pts[static_cast<std::size_t>(i)].x,
                                     centers[static_cast<std::size_t>(j)], sigma);
        y(i) = pts[static_cast<std::size_t>(i)].y;
    }
    const Eigen::VectorXd c = A.colPivHouseholderQr().solve(y);

    coefs.reserve(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) coefs.push_back(c(j));
    return coefs;
}

double evaluateGaussLeastSquares(const std::vector<double>& coefs,
                                 const std::vector<double>& centers,
                                 double sigma, double x) {
    const int m = static_cast<int>(centers.size());
    if (coefs.size() != centers.size()) return 0.0;
    double y = 0.0;
    for (int j = 0; j < m; ++j)
        y += coefs[static_cast<std::size_t>(j)] *
             gaussNodeValue(x, centers[static_cast<std::size_t>(j)], sigma);
    return y;
}

std::vector<Point> gaussLeastSquares(const std::vector<Point>& pts, double sigma,
                                     int numCenters, int samples) {
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    std::vector<double> centers;
    const std::vector<double> coefs = fitGaussLeastSquares(pts, sigma, numCenters, centers);
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
        out.push_back({x, evaluateGaussLeastSquares(coefs, centers, sigma, x)});
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

// ---------------------------------------------------------------------------
// 单参数曲线拟合（作业3）
// ---------------------------------------------------------------------------

std::vector<double> parameterizeUniform(const std::vector<Point>& pts) {
    std::vector<double> t(pts.size());
    double inv_size = 1.0 / pts.size();
    for (std::size_t i = 0; i < pts.size(); ++i) t[i] = static_cast<double>(i) * inv_size;
    return t;
}

std::vector<double> parameterizeChordal(const std::vector<Point>& pts)
{
    std::vector<double> t(pts.size());
    double chord_sum = 0;
    t[0] = 0;
    for (std::size_t i = 1; i < pts.size(); ++i)
    {
        chord_sum += L2(pts[i - 1], pts[i]);
        t[i] = chord_sum;
    }
    double inv_chord_sum = 1.0 / chord_sum;
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        t[i] *= inv_chord_sum;
    }

    return t;
}

std::vector<double> parameterizeCentripetal(const std::vector<Point>& pts)
{
    std::vector<double> t(pts.size());
    double len_sum = 0;
    t[0] = 0;
    for (std::size_t i = 1; i < pts.size(); ++i)
    {
        len_sum += std::sqrt(L2(pts[i - 1], pts[i]));
        t[i] = len_sum;
    }
    double inv_len_sum = 1.0 / len_sum;
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        t[i] *= inv_len_sum;
    }

    return t;
}

std::vector<double> parameterizeFoley(const std::vector<Point>& pts)
{
    std::vector<double> t(pts.size());
    if (pts.size() < 2) return t;  // 防御：空/单点直接返回全 0

    double len_sum = 0;
    t[0] = 0;
    // t1 use pts0, pts0, pts1 to calculate len

    for (std::size_t i = 1; i < pts.size(); ++i)
    {
        // 注意：i 是无符号类型，直接写 i - 2 在 i < 2 时会下溢成巨大数，
        // 导致 pts[巨大数] 越界；必须先判断再减。
        const std::size_t i0 = (i >= 2) ? (i - 2) : 0;
        Point p0 = pts[i0];
        Point p1 = pts[i - 1];
        Point p2 = pts[i]; // current t 
        Point p3 = pts[std::min(i + 1, pts.size() - 1)];
        Eigen::Vector2d v1(p1.x - p0.x, p1.y - p0.y);
        Eigen::Vector2d v2(p2.x - p1.x, p2.y - p1.y);
        Eigen::Vector2d v3(p3.x - p2.x, p3.y - p2.y);

        double p0p1L2 = L2(p0, p1);
        double p1p2L2 = L2(p1, p2);
        double p2p3L2 = L2(p2, p3);

        double alpha_0 = std::min(EIGEN_PI - angleInRadians(v1, v2), EIGEN_PI / 2);
        double alpha_1 = std::min(EIGEN_PI - angleInRadians(v2, v3), EIGEN_PI / 2);

        double len = p1p2L2 * (1 + 1.5 * alpha_0 * p0p1L2 / (p0p1L2 + p1p2L2) + 1.5 * alpha_1 * p1p2L2 / (p1p2L2 + p2p3L2));
        len_sum += len;
        t[i] = len_sum;
    }
    if (!(len_sum > 0.0)) {
        // 所有点重合（弦长全 0）或计算中出现 NaN：返回干净的全 0，
        // 不能直接 return t（此时 t 里可能已残留 NaN）。
        std::fill(t.begin(), t.end(), 0.0);
        return t;
    }
    double inv_len_sum = 1.0 / len_sum;
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        t[i] *= inv_len_sum;
    }

    return t;
}

std::vector<Point> fitParametricCurve(const std::vector<Point>& pts,
                                      const std::vector<double>& t,
                                      int degree,
                                      int samples) {
    std::vector<Point> out;
    if (pts.size() < 2 || t.size() != pts.size()) return out;

    // 1) 构造 (t, x) 和 (t, y) 两个标量点列 —— 作业要求的"构造新的点列"
    std::vector<Point> px, py;
    px.reserve(pts.size());
    py.reserve(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        px.push_back({t[i], pts[i].x});
        py.push_back({t[i], pts[i].y});
    }

    // 2) 分别做 degree 次最小二乘拟合，得到 x(t)、y(t) 的采样点列。
    //    两次采样用的 t 网格相同（都取自各自数据点的 t_min..t_max），
    //    因此可以按下标直接合成二维点。
    const std::vector<Point> cx = leastSquaresPolynomial(px, degree, samples);
    const std::vector<Point> cy = leastSquaresPolynomial(py, degree, samples);
    if (cx.empty() || cy.empty()) return out;

    // 3) 合成二维曲线 (x(t_k), y(t_k))
    out.reserve(cx.size());
    for (std::size_t k = 0; k < cx.size(); ++k) out.push_back({cx[k].y, cy[k].y});
    return out;
}

std::vector<Point> fitParametricCurveGauss(const std::vector<Point>& pts,
                                           const std::vector<double>& t,
                                           double sigma,
                                           int numCenters,
                                           int samples) {
    std::vector<Point> out;
    if (pts.size() < 2 || t.size() != pts.size()) return out;

    // 与幂基版完全同构：拆 (t,x) (t,y) 两个标量点列，只是分别改用 Gauss 基最小二乘
    std::vector<Point> px, py;
    px.reserve(pts.size());
    py.reserve(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        px.push_back({t[i], pts[i].x});
        py.push_back({t[i], pts[i].y});
    }

    const std::vector<Point> cx = gaussLeastSquares(px, sigma, numCenters, samples);
    const std::vector<Point> cy = gaussLeastSquares(py, sigma, numCenters, samples);
    if (cx.empty() || cy.empty()) return out;

    out.reserve(cx.size());
    for (std::size_t k = 0; k < cx.size(); ++k) out.push_back({cx[k].y, cy[k].y});
    return out;
}

// 内部：标量三次样条采样器（无切线控制版本，保留作兼容/测试用；要求 pts 的 x 严格递增）
std::vector<Point> cubicSpline(const std::vector<Point>& pts, CubicSplineType splineType,
                               int perSegmentSamples)
{
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    const std::vector<double> coefs = cubicSplineCoefs(pts, splineType);
    if (coefs.empty()) return out;

    perSegmentSamples = std::max(perSegmentSamples, 2);
    const int segments = static_cast<int>(pts.size()) - 1;
    out.reserve(static_cast<std::size_t>(perSegmentSamples * segments));

    for (int i = 0; i < segments; i++)
    {
        double xmin = pts[i].x < pts[i + 1].x ? pts[i].x : pts[i + 1].x;
        double xmax = pts[i].x < pts[i + 1].x ? pts[i + 1].x : pts[i].x;
        std::vector<double> segmentCoefs = { coefs[4 * i], coefs[4 * i + 1], coefs[4 * i + 2], coefs[4 * i + 3] };
        for (int k = 0; k < perSegmentSamples; ++k)
        {
            const double x = xmin + (xmax - xmin) * static_cast<double>(k) / (perSegmentSamples - 1);
            out.push_back({ x, evaluateCubicSplineSegment(pts[i], segmentCoefs, x) });
        }
    }

    return out;
}

namespace {
// 每节点的一阶导（标量系统）：自由节点 mLeft == mRight == 未知数 m_j；
// 被编辑节点（钉值）mLeft/mRight 是给定数据（左右可不同，最通用）。
struct NodeSlope {
    bool free = true;
    double mLeft = 0.0;
    double mRight = 0.0;
};
}  // namespace

// 三转角法标量引擎（Hermite 记法）：每节点存左右一阶导 (mLeft, mRight)，
// 自由节点左右相等（一个未知数 m_j）。方程对每个自由节点 j：
//   内部 C²：λ_j·m_{j−1}^R + 2λ_j·m_j^L + 2μ_j·m_j^R + μ_j·m_{j+1}^L
//            = 3[ λ_j·Δy_{j−1}/h_{j−1} + μ_j·Δy_j/h_j ]
//    其中 λ_j = h_j/(h_{j−1}+h_j)、μ_j = h_{j−1}/(h_{j−1}+h_j)；
//    全自由时退化为教科书三转角方程 λ_j·m_{j−1} + 2m_j + μ_j·m_{j+1} = g_j。
//   自然端（n+1 个方程）：2m₀^R + m₁^L = 3Δy₀/h₀；mₙ₋₁^R + 2mₙ^L = 3Δyₙ₋₁/hₙ₋₁
//   夹持端（n−1 个方程）：m₀、mₙ 已知（端点斜率相邻差分估算），只解内部节点。
// pins 为空或长度不符 = 全部自由；pins[k].free=false 时左右切线为钉值。
// 返回每节点最终 (mLeft, mRight)。
std::vector<NodeSlope> cubicSplineNodeSlopes(const std::vector<Point>& pts,
                                             CubicSplineType splineType,
                                             const std::vector<NodeSlope>& pins)
{
    const int n = static_cast<int>(pts.size()) - 1;
    if (n < 1) return {};
    for (int i = 0; i < n; ++i)
        if (pts[i + 1].x <= pts[i].x) return {};

    const bool hasPins = pins.size() == pts.size();

    // 是否自由：夹持端点 m₀、mₙ 一律视为钉值（得到 n−1 个内部自由节点）
    const auto isFree = [&](int k) {
        if (hasPins && !pins[static_cast<std::size_t>(k)].free) return false;
        if (splineType == CubicSplineType::Clamped && (k == 0 || k == n)) return false;
        return true;
    };
    // 钉值：编辑值优先，否则夹持端点用相邻差分估算（仅对非自由节点调用）
    const auto pinVal = [&](int k, bool rightSide) -> double {
        if (hasPins && !pins[static_cast<std::size_t>(k)].free)
            return rightSide ? pins[static_cast<std::size_t>(k)].mRight
                             : pins[static_cast<std::size_t>(k)].mLeft;
        if (k == 0) return (pts[1].y - pts[0].y) / (pts[1].x - pts[0].x);
        if (k == n) return (pts[n].y - pts[n - 1].y) / (pts[n].x - pts[n - 1].x);
        return 0.0;
    };

    // 自由未知量映射
    std::vector<int> nodeToFree(pts.size(), -1);
    std::vector<int> freeNodes;
    for (int k = 0; k <= n; ++k)
        if (isFree(k)) {
            nodeToFree[static_cast<std::size_t>(k)] = static_cast<int>(freeNodes.size());
            freeNodes.push_back(k);
        }
    const int numFree = static_cast<int>(freeNodes.size());

    std::vector<double> mFree(static_cast<std::size_t>(numFree), 0.0);
    if (numFree > 0) {
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(numFree, numFree);
        Eigen::VectorXd y = Eigen::VectorXd::Zero(numFree);

        // (coeff × 节点 k 某侧切线) 入系统：自由 → A 列；钉值 → 右端项
        const auto addTerm = [&](int row, int k, bool rightSide, double coeff) {
            if (isFree(k)) {
                A(row, nodeToFree[static_cast<std::size_t>(k)]) += coeff;
            } else {
                y(row) -= coeff * pinVal(k, rightSide);
            }
        };

        for (int row = 0; row < numFree; ++row) {
            const int j = freeNodes[static_cast<std::size_t>(row)];
            if (j == 0) {
                // 自然端 x₀：2m₀^R + m₁^L = 3Δy₀/h₀
                const double h0 = pts[1].x - pts[0].x;
                const double dy0 = pts[1].y - pts[0].y;
                addTerm(row, 0, true, 2.0);
                addTerm(row, 1, false, 1.0);
                y(row) += 3.0 * dy0 / h0;
            } else if (j == n) {
                // 自然端 xₙ：mₙ₋₁^R + 2mₙ^L = 3Δyₙ₋₁/hₙ₋₁
                const double hn = pts[n].x - pts[n - 1].x;
                const double dyn = pts[n].y - pts[n - 1].y;
                addTerm(row, n - 1, true, 1.0);
                addTerm(row, n, false, 2.0);
                y(row) += 3.0 * dyn / hn;
            } else {
                // 内部 C²（Hermite 记法，自由节点左右相等时退化为 λ_j m_{j−1} + 2m_j + μ_j m_{j+1} = g_j）
                const double hPrev = pts[j].x - pts[j - 1].x;   // h_{j-1}
                const double hNext = pts[j + 1].x - pts[j].x;   // h_j
                const double dyPrev = pts[j].y - pts[j - 1].y;  // Δy_{j-1}
                const double dyNext = pts[j + 1].y - pts[j].y;  // Δy_j
                const double sum = hPrev + hNext;
                const double lam = hNext / sum;  // λ_j
                const double mu = hPrev / sum;   // μ_j
                addTerm(row, j - 1, true, lam);
                addTerm(row, j, false, 2.0 * lam);
                addTerm(row, j, true, 2.0 * mu);
                addTerm(row, j + 1, false, mu);
                y(row) += 3.0 * (lam * dyPrev / hPrev + mu * dyNext / hNext);
            }
        }

        const Eigen::VectorXd sol = A.lu().solve(y);
        for (int j = 0; j < numFree; ++j)
            mFree[static_cast<std::size_t>(j)] = sol(j);
    }

    // 组装每节点的左右切线
    std::vector<NodeSlope> out(pts.size());
    for (int k = 0; k <= n; ++k) {
        if (isFree(k)) {
            out[static_cast<std::size_t>(k)].free = true;
            out[static_cast<std::size_t>(k)].mLeft = out[static_cast<std::size_t>(k)].mRight =
                mFree[static_cast<std::size_t>(nodeToFree[static_cast<std::size_t>(k)])];
        } else {
            out[static_cast<std::size_t>(k)].free = false;
            out[static_cast<std::size_t>(k)].mLeft = pinVal(k, false);
            out[static_cast<std::size_t>(k)].mRight = pinVal(k, true);
        }
    }
    return out;
}

// Hermite 分段插值求值（三转角法 / 三次样条的段表示基础）：
// 第 j 段 [x_j, x_{j+1}] 长 h，左端 (f_j, m_j^R)、右端 (f_{j+1}, m_{j+1}^L)，
// u=(x−x_j)/h：S(u) = H00·f_j + H10·(h·m_j^R) + H01·f_{j+1} + H11·(h·m_{j+1}^L)
double hermiteEvaluate(double x, double xj, double h, double fj, double mjR,
                       double fj1, double mj1L)
{
    if (h <= 0.0) return 0.0;
    const double u = (x - xj) / h;
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double H00 = 2 * u3 - 3 * u2 + 1;  // (1+2u)(1−u)²
    const double H10 = u3 - 2 * u2 + u;      // u(u−1)²
    const double H01 = -2 * u3 + 3 * u2;     // u²(3−2u)
    const double H11 = u3 - u2;              // u²(u−1)
    return H00 * fj + H10 * (h * mjR) + H01 * fj1 + H11 * (h * mj1L);
}

// 内部：按 Hermite 逐段采样（段 j 用 m_j^R、m_{j+1}^L，最通用的表示）
std::vector<Point> sampleHermiteSpline(const std::vector<Point>& pts,
                                       const std::vector<NodeSlope>& nodes,
                                       int perSegmentSamples)
{
    std::vector<Point> out;
    const int n = static_cast<int>(pts.size()) - 1;
    if (n < 1 || nodes.size() != pts.size()) return out;
    perSegmentSamples = std::max(perSegmentSamples, 2);
    out.reserve(static_cast<std::size_t>(perSegmentSamples * n));
    for (int j = 0; j < n; ++j) {
        const double xj = pts[j].x;
        const double h = pts[j + 1].x - xj;
        for (int k = 0; k < perSegmentSamples; ++k) {
            const double x = xj + h * static_cast<double>(k) / (perSegmentSamples - 1);
            out.push_back({x, hermiteEvaluate(x, xj, h, pts[j].y, nodes[j].mRight,
                                              pts[j + 1].y, nodes[j + 1].mLeft)});
        }
    }
    return out;
}

// 兼容层：由节点切线换算成每段幂基系数 [a,b,c,d]
//（与 Hermite 形式数学等价：b=m_j^R，c=(3Δy−h(2m_j^R+m_{j+1}^L))/h²，
//  d=(−2Δy+h(m_j^R+m_{j+1}^L))/h³）
std::vector<double> cubicSplineCoefs(const std::vector<Point>& pts, CubicSplineType splineType)
{
    const std::vector<NodeSlope> nodes = cubicSplineNodeSlopes(pts, splineType, {});
    if (nodes.size() != pts.size()) return {};
    const int n = static_cast<int>(pts.size()) - 1;
    std::vector<double> coefs;
    coefs.reserve(static_cast<std::size_t>(4 * n));
    for (int j = 0; j < n; ++j) {
        const double h = pts[j + 1].x - pts[j].x;
        const double dy = pts[j + 1].y - pts[j].y;
        const double l = nodes[j].mRight;
        const double r = nodes[j + 1].mLeft;
        coefs.push_back(pts[j].y);
        coefs.push_back(l);
        coefs.push_back((3 * dy - h * (2 * l + r)) / (h * h));
        coefs.push_back((-2 * dy + h * (l + r)) / (h * h * h));
    }
    return coefs;
}
// this only accepts 4 coefs
double evaluateCubicSplineSegment(const Point& stP, const std::vector<double>& coefs, double x)
{
    if (coefs.size() != 4) return 0;
    double rst = 0;
    double power = 1;
    
    for (int i = 0; i < 4; i++)
    {
        rst += coefs[i] * power;
        power *= (x - stP.x);
    }

    return rst;
}

// 参数型三次样条（三转角法 + 弦长参数化）——三次样条的通用入口。
// 对任意有序点列先用弦长参数化得到 tᵢ（严格递增），再对 (tᵢ,xᵢ)、(tᵢ,yᵢ)
// 分别做标量三次样条并合成二维曲线，因此 x、y 无需单调。
// controls（可选，长度须等于 pts 数量）：per-node 切线覆盖，见 VertexTangent；
// 被编辑节点（Smooth/Straight/Corner）的切线作为数据（左右一阶导），
// 该节点不再要求 C²。
std::vector<Point> fitParametricSpline(const std::vector<Point>& pts,
                                       CubicSplineType splineType,
                                       int perSegmentSamples,
                                       const std::vector<VertexTangent>& controls)
{
    std::vector<Point> out;
    if (pts.size() < 2) return out;

    const std::vector<double> t = parameterizeChordal(pts);
    if (t.size() != pts.size()) return out;

    // 1) 拆 (t, x) 和 (t, y) 两个标量点列
    std::vector<Point> px, py;
    px.reserve(pts.size());
    py.reserve(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        px.push_back({t[i], pts[i].x});
        py.push_back({t[i], pts[i].y});
    }

    // 2) 2D 切线控制 → x/y 两个标量系统的每节点 (free, mLeft, mRight)
    //    （Smooth 左右相等；Straight/Corner 左右独立，都是钉值）
    const bool hasCtrl = controls.size() == pts.size();
    const auto makePins = [&](bool forX) {
        std::vector<NodeSlope> pins;
        if (!hasCtrl) return pins;
        pins.resize(pts.size());
        for (std::size_t k = 0; k < pts.size(); ++k) {
            const VertexTangent& c = controls[k];
            NodeSlope& p = pins[k];
            if (c.mode == VertexTangent::Mode::Free) {
                p.free = true;
            } else if (c.mode == VertexTangent::Mode::Smooth) {
                p.free = false;
                p.mLeft = p.mRight = forX ? c.t.x : c.t.y;
            } else {  // Straight / Corner
                p.free = false;
                p.mLeft = forX ? c.tLeft.x : c.tLeft.y;
                p.mRight = forX ? c.tRight.x : c.tRight.y;
            }
        }
        return pins;
    };
    const std::vector<NodeSlope> pinsX = makePins(true);
    const std::vector<NodeSlope> pinsY = makePins(false);

    // 3) 分别求各节点左右切线，并按 Hermite 逐段采样；
    //    两边 t 网格相同（同参数化、同每段采样数），可按下标合成。
    const std::vector<NodeSlope> nx = cubicSplineNodeSlopes(px, splineType, pinsX);
    const std::vector<NodeSlope> ny = cubicSplineNodeSlopes(py, splineType, pinsY);
    if (nx.size() != pts.size() || ny.size() != pts.size()) return out;

    const std::vector<Point> sx = sampleHermiteSpline(px, nx, perSegmentSamples);
    const std::vector<Point> sy = sampleHermiteSpline(py, ny, perSegmentSamples);
    if (sx.empty() || sy.empty()) return out;

    // 4) 合成二维曲线 (x(t_k), y(t_k))
    out.reserve(sx.size());
    for (std::size_t k = 0; k < sx.size(); ++k) out.push_back({sx[k].y, sy[k].y});
    return out;
}

// ---------------------------------------------------------------------------
// 分段三次 Bezier（Catmull-Rom 控制点构造）
// ---------------------------------------------------------------------------

double binomialCoef(unsigned n, unsigned i)
{
    long long numerator = 1, denominator = 1;
    for (unsigned k = n; k > n - i; k--) numerator *= k;
    for (unsigned k = i; k > 0; k--) denominator *= k;
    return double(numerator) / denominator;
}

// 三次 Bernstein 基求值：B(t) = Σ Bᵢ·C(3,i)·t^i·(1−t)^(3−i)，pts 恰为 4 个控制点。
// 注意 (1−t) 的指数是 3−i，与二项式系数 C(3,i) 匹配（早期版本误写为 4−i，
// 导致幂次错位、曲线形状错误，已修正）。
Point evaluateCubicBernsteinBasisCurve(const std::vector<Point>& pts, double t)
{
    if (pts.size() != 4) return {0, 0};

    double rstX = 0;
    double rstY = 0;
    for (int i = 0; i < 4; i++)
    {
        const double w = binomialCoef(3, i) * std::pow(t, i) * std::pow(1 - t, 3 - i);
        rstX += pts[i].x * w;
        rstY += pts[i].y * w;
    }
    return { rstX, rstY };
}

// Catmull-Rom 虚拟点：index 越界（-1 / n）时用端点镜像外推补齐邻居，
// 使首尾段也能套用统一的切线公式；范围内则原样返回数据点。
Point catmullRomVirtualPoint(const std::vector<Point>& pts, int index)
{
    const int n = static_cast<int>(pts.size());
    if (n < 2) return pts.empty() ? Point{0, 0} : pts[0];
    if (index < 0) {
        // P_{−1} = 2·P₀ − P₁（把 P₁ 关于 P₀ 镜像到 P₀ 另一侧）
        return { 2 * pts[0].x - pts[1].x, 2 * pts[0].y - pts[1].y };
    }
    if (index >= n) {
        // P_n = 2·P_{n−1} − P_{n−2}（把 P_{n−2} 关于 P_{n−1} 镜像）
        return { 2 * pts[n - 1].x - pts[n - 2].x, 2 * pts[n - 1].y - pts[n - 2].y };
    }
    return pts[index];
}

BezierSegment buildBezierSegment(const std::vector<Point>& pts, int i)
{
    BezierSegment seg;
    const int n = static_cast<int>(pts.size());
    if (n < 2 || i < 0 || i >= n - 1) return seg;  // 非法段：返回空结构（全 0）

    if (n == 2) {
        // 只有一段、没有相邻点：不需要创建虚拟点，B1/B2 直接取两端型值点 → 直线段
        seg.control[0] = pts[0];
        seg.control[1] = pts[1];
        return seg;
    }

    // Catmull-Rom 切线 → Bezier 中间控制点（边界段用虚拟点补齐邻居）：
    //   B1 = P_i + (P_{i+1} − P_{i−1})/6，B2 = P_{i+1} − (P_{i+2} − P_i)/6
    const Point pm1 = (i == 0) ? catmullRomVirtualPoint(pts, -1) : pts[i - 1];
    const Point pp2 = (i == n - 2) ? catmullRomVirtualPoint(pts, n) : pts[i + 2];
    seg.control[0] = { pts[i].x + (pts[i + 1].x - pm1.x) / 6.0,
                       pts[i].y + (pts[i + 1].y - pm1.y) / 6.0 };
    seg.control[1] = { pts[i + 1].x - (pp2.x - pts[i].x) / 6.0,
                       pts[i + 1].y - (pp2.y - pts[i].y) / 6.0 };
    return seg;
}

std::vector<BezierSegment> buildBezierSegments(const std::vector<Point>& pts)
{
    std::vector<BezierSegment> segs;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return segs;
    segs.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 0; i < n - 1; ++i)
        segs.push_back(buildBezierSegment(pts, i));
    return segs;
}

// 第 i 段（P_i → P_{i+1}）在参数 u 处求值：
// 结合本段 segment 的两个中间控制点与自身两个型值点 pts[i]、pts[i+1] 做三次插值。
Point evaluateBezierSegment(const std::vector<Point>& pts, int i,
                            const BezierSegment& seg, double u)
{
    const int n = static_cast<int>(pts.size());
    if (n < 2 || i < 0 || i >= n - 1) return {0, 0};
    u = std::clamp(u, 0.0, 1.0);
    // 4 个控制点 = 本段中间控制点 + 两端型值点
    const Point ctrl[4] = { pts[i], seg.control[0], seg.control[1], pts[i + 1] };
    return evaluateCubicBernsteinBasisCurve(std::vector<Point>(ctrl, ctrl + 4), u);
}

std::vector<Point> bezierCatmullRomInterpolate(const std::vector<Point>& pts,
                                               int perSegmentSamples)
{
    // 不带段控制点：segments 为空 → 内部自动按 Catmull-Rom 重建
    return bezierCatmullRomInterpolate(pts, {}, perSegmentSamples);
}

std::vector<Point> bezierCatmullRomInterpolate(const std::vector<Point>& pts,
                                               const std::vector<BezierSegment>& segs,
                                               int perSegmentSamples)
{
    std::vector<Point> out;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return out;

    // 调用方给的段控制点数量不对（如 GUI 拖过之后 pts 变了）→ 回退自动重建
    const std::vector<BezierSegment>& effective =
        (segs.size() == static_cast<std::size_t>(n - 1)) ? segs : buildBezierSegments(pts);

    perSegmentSamples = std::max(perSegmentSamples, 2);
    out.reserve(static_cast<std::size_t>(perSegmentSamples * (n - 1)));
    for (int i = 0; i < n - 1; ++i) {
        for (int k = 0; k < perSegmentSamples; ++k) {
            const double u = static_cast<double>(k) / (perSegmentSamples - 1);
            out.push_back(evaluateBezierSegment(pts, i, effective[static_cast<std::size_t>(i)], u));
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Chaikin 细分（割角法，逼近型）
// ---------------------------------------------------------------------------

// 单轮细分：每条边 (a, b) 生成两个割角点 3/4·a + 1/4·b 与 1/4·a + 3/4·b，
// 老点全部抛弃；开放曲线保留首尾端点（端点不割角）。
std::vector<Point> chaikinSubdivideOnce(const std::vector<Point>& pts)
{
    std::vector<Point> out;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return out;
    if (n == 2) return pts;  // 两个点本身就是直线段，无需割角

    out.reserve(static_cast<std::size_t>(2 * n));
    out.push_back(pts.front());  // 保留起点
    for (int i = 0; i + 1 < n; ++i) {
        const Point& a = pts[i];
        const Point& b = pts[i + 1];
        out.push_back({0.75 * a.x + 0.25 * b.x, 0.75 * a.y + 0.25 * b.y});  // ν′_{2i+1}
        out.push_back({0.25 * a.x + 0.75 * b.x, 0.25 * a.y + 0.75 * b.y});  // ν′_{2i+2}
    }
    out.push_back(pts.back());  // 保留终点
    return out;
}

// 细分曲线：对 pts 迭代 iterations 轮，把最终多边形直接作为逼近曲线返回。
//（迭代骨架复用 subdivideIterate，与均匀三次 B 样条细分共用。）
std::vector<Point> chaikinSubdivisionCurve(const std::vector<Point>& pts, int iterations)
{
    return subdivideIterate(pts, iterations, &chaikinSubdivideOnce);
}

// ---------------------------------------------------------------------------
// 均匀三次 B 样条细分（逼近型）
// ---------------------------------------------------------------------------

// 单轮细分：内部顶点平滑 ν′_{2i} = 1/8·ν_{i−1} + 3/4·ν_i + 1/8·ν_{i+1}，
// 每条边插中点 ν′_{2i+1} = 1/2·ν_i + 1/2·ν_{i+1}；老点抛弃，首尾端点保留（钳制端）。
std::vector<Point> bsplineSubdivideOnce(const std::vector<Point>& pts)
{
    std::vector<Point> out;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return out;
    if (n == 2) return pts;  // 两个点本身就是直线段，无需细分

    out.reserve(static_cast<std::size_t>(2 * n - 1));
    out.push_back(pts.front());  // 保留起点 ν′₀ = ν₀
    for (int i = 0; i + 1 < n; ++i) {
        const Point& a = pts[i];
        const Point& b = pts[i + 1];
        // 边 (a, b) 的中点 ν′_{2i+1}
        out.push_back({0.5 * a.x + 0.5 * b.x, 0.5 * a.y + 0.5 * b.y});
        // 若 b 是内部顶点（非最后一个），做 1/8:3/4:1/8 平滑 ν′_{2i+2}
        if (i + 1 < n - 1) {
            const Point& c = pts[i + 2];
            out.push_back({0.125 * a.x + 0.75 * b.x + 0.125 * c.x,
                           0.125 * a.y + 0.75 * b.y + 0.125 * c.y});
        }
    }
    out.push_back(pts.back());  // 保留终点 ν′_{2n−1} = ν_{n−1}
    return out;
}

// 细分曲线：迭代骨架与 Chaikin 完全共用（subdivideIterate），只换单轮规则。
std::vector<Point> bsplineSubdivisionCurve(const std::vector<Point>& pts, int iterations)
{
    return subdivideIterate(pts, iterations, &bsplineSubdivideOnce);
}

std::vector<Point> fourPointSubdivideOnce(const std::vector<Point>& pts)
{
    std::vector<Point> out;
    const int n = static_cast<int>(pts.size());
    if (n < 2) return out;

    // 2个点：直接返回（或可改为插入中点，但这里按你的逻辑返回原样）
    if (n == 2) return pts;

    out.reserve(static_cast<std::size_t>(2 * n - 1));

    // 保留第一个旧点
    out.push_back(pts[0]);

    for (int i = 0; i + 1 < n; ++i)
    {
        // 获取4个点：v_{i-1}, v_i, v_{i+1}, v_{i+2}
        // 使用镜像点处理边界

        Point v_im1;  // v_{i-1}
        if (i == 0)
        {
            // 镜像：v_{-1} = 2*v_0 - v_1
            v_im1.x = 2.0 * pts[0].x - pts[1].x;
            v_im1.y = 2.0 * pts[0].y - pts[1].y;
        }
        else
        {
            v_im1 = pts[i - 1];
        }

        const Point& v_i = pts[i];
        const Point& v_ip1 = pts[i + 1];

        Point v_ip2;  // v_{i+2}
        if (i + 2 >= n)
        {
            // 镜像：v_{n} = 2*v_{n-1} - v_{n-2}
            v_ip2.x = 2.0 * pts[n - 1].x - pts[n - 2].x;
            v_ip2.y = 2.0 * pts[n - 1].y - pts[n - 2].y;
        }
        else
        {
            v_ip2 = pts[i + 2];
        }

        // 4点插值公式计算新插入的点
        // v'_{2i+1} = 9/16*(v_i + v_{i+1}) - 1/16*(v_{i-1} + v_{i+2})
        Point newPoint;
        newPoint.x = (9.0 / 16.0) * (v_i.x + v_ip1.x) - (1.0 / 16.0) * (v_im1.x + v_ip2.x);
        newPoint.y = (9.0 / 16.0) * (v_i.y + v_ip1.y) - (1.0 / 16.0) * (v_im1.y + v_ip2.y);

        out.push_back(newPoint);  // 插入新点

        // 保留旧点 v_{i+1}
        out.push_back(v_ip1);
    }

    // 注意：由于循环中每次都会push_back(v_ip1)，最后一个旧点已经被保留
    // 所以不需要再额外push_back(pts.back())
    // 但如果n==2，已经在前面返回了

    return out;  // 长度为 2n-1
}

std::vector<Point> fourPointSubdivisionCurve(const std::vector<Point>& pts, int iterations)
{
    return subdivideIterate(pts, iterations, &fourPointSubdivideOnce);
}

}  // namespace curve
