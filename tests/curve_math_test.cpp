// 数学库自测：验证曲线算法的数值正确性。
// 纯 C++、无 Qt 依赖，构建后直接运行 exe 即可：
//   cmake --build build\msvc2022 --config Debug --target curve_math_test
//   build\msvc2022\Debug\curve_math_test.exe
#include "math/CurveMath.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {
int g_failures = 0;

void check(bool cond, const char* expr, int line) {
    if (!cond) {
        std::printf("FAIL  line %d: %s\n", line, expr);
        ++g_failures;
    }
}
}  // namespace

#define CHECK(cond) check((cond), #cond, __LINE__)

int main() {
    // ---------- 1) Lagrange 插值必须经过每个数据点 ----------
    {
        const std::vector<curve::Point> pts = {{0.1, 0.2}, {0.3, 0.9}, {0.7, 0.4}, {0.95, 0.6}};
        for (const auto& p : pts)
            CHECK(std::abs(curve::lagrangeAt(pts, p.x) - p.y) < 1e-9);

        const auto interp = curve::lagrangeInterpolate(pts, 100);
        CHECK(interp.size() == 100);
        // 采样点两端的 x 必须落在数据范围内（浮点运算，用容差比较）
        CHECK(std::abs(interp.front().x - pts.front().x) < 1e-12);
        CHECK(std::abs(interp.back().x - pts.back().x) < 1e-12);
    }

    // ---------- 2) 最小二乘：能精确恢复抛物线 y = 2 + 3x - 4x^2 ----------
    {
        std::vector<curve::Point> para;
        for (int i = 0; i <= 10; ++i) {
            const double x = 0.1 * i;
            para.push_back({x, 2.0 + 3.0 * x - 4.0 * x * x});
        }
        const auto c = curve::fitPolynomial(para, 2);
        CHECK(c.size() == 3);
        CHECK(std::abs(c[0] - 2.0) < 1e-8);
        CHECK(std::abs(c[1] - 3.0) < 1e-8);
        CHECK(std::abs(c[2] + 4.0) < 1e-8);

        const auto approx = curve::leastSquaresPolynomial(para, 2, 50);
        CHECK(approx.size() == 50);
        for (const auto& p : approx)
            CHECK(std::abs(curve::evaluatePolynomial(c, p.x) - p.y) < 1e-8);
    }

    // ---------- 3) 次数自动限制：degree > n-1 时按 n-1 处理 ----------
    {
        std::vector<curve::Point> pts;
        for (int i = 0; i <= 4; ++i) pts.push_back({0.2 * i, std::sin(0.2 * i)});
        const auto c = curve::fitPolynomial(pts, 100);  // n=5 -> 最多 4 次
        CHECK(c.size() == 5);
    }

    // ---------- 4) 过 3 点直线：二次拟合应精确经过每个点 ----------
    {
        const std::vector<curve::Point> tri = {{0.0, 1.0}, {0.5, 2.0}, {1.0, 3.0}};
        const auto ct = curve::fitPolynomial(tri, 2);
        CHECK(ct.size() == 3);
        for (const auto& p : tri)
            CHECK(std::abs(curve::evaluatePolynomial(ct, p.x) - p.y) < 1e-9);
    }

    // ---------- 5) 重复 x 不崩溃（Eigen QR 的鲁棒性） ----------
    {
        const std::vector<curve::Point> dup = {{0.2, 0.3}, {0.2, 0.9}, {0.6, 0.5}};
        const auto c = curve::fitPolynomial(dup, 2);
        CHECK(!c.empty());
    }

    // ---------- 6) 幂基插值：精确恢复原多项式，且与 Lagrange 结果一致 ----------
    {
        // 4 个点取自抛物线 y = 2 + 3x - 4x^2
        const auto poly = [](double x) { return 2.0 + 3.0 * x - 4.0 * x * x; };
        const std::vector<curve::Point> pts = {
            {0.1, poly(0.1)}, {0.3, poly(0.3)}, {0.7, poly(0.7)}, {0.95, poly(0.95)}};

        const auto pb = curve::powerBasePolynomialInterpolate(pts, 200);
        CHECK(pb.size() == 200);
        // 插值多项式唯一，且数据来自该抛物线，所以每个采样点都应落在抛物线上
        for (const auto& p : pb) CHECK(std::abs(p.y - poly(p.x)) < 1e-8);

        // 与 Lagrange 插值画出的曲线一致（同一多项式、不同表示）
        const auto lg = curve::lagrangeInterpolate(pts, 200);
        CHECK(lg.size() == pb.size());
        for (std::size_t i = 0; i < pb.size(); ++i)
            CHECK(std::abs(lg[i].y - pb[i].y) < 1e-8);
    }

    // ---------- 7) 高斯基函数插值：必须精确经过每一个数据点（含最后一个） ----------
    {
        const double sigma = 0.3;
        std::vector<curve::Point> pts;
        for (int i = 0; i < 5; ++i) {
            const double x = 0.15 + 0.18 * i;  // 0.15 ~ 0.87
            pts.push_back({x, std::sin(3.0 * x) + 0.2 * x});
        }

        const auto coefs = curve::fitGaussBasePolynomial(pts, sigma);
        CHECK(coefs.size() == pts.size() + 1);  // n 个高斯系数 + 1 个常数项
        // 在每个数据点上求值都应精确等于 y（重点：最后一个点 i = n-1）
        for (std::size_t i = 0; i < pts.size(); ++i)
            CHECK(std::abs(curve::evaluateGaussBasePolynomial(coefs, pts, sigma, pts[i].x) -
                           pts[i].y) < 1e-6);

        const auto gaussCurve = curve::gaussBasePolynomialInterpolate(pts, sigma, 200);
        CHECK(gaussCurve.size() == 200);
        CHECK(std::abs(gaussCurve.front().x - pts.front().x) < 1e-12);
        CHECK(std::abs(gaussCurve.back().x - pts.back().x) < 1e-12);
    }

    if (g_failures == 0) {
        std::printf("ALL TESTS PASSED\n");
        return 0;
    }
    std::printf("%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
