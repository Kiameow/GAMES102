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

    // ---------- 8) 均匀参数化：tᵢ = i/n（归一化到 (0,1]），严格递增 ----------
    {
        const std::vector<curve::Point> pts = {
            {0.0, 0.0}, {0.2, 1.0}, {0.8, 0.5}, {1.5, 2.0}};
        const auto t = curve::parameterizeUniform(pts);
        CHECK(t.size() == pts.size());
        for (std::size_t i = 0; i < t.size(); ++i) {
            CHECK(t[i] == static_cast<double>(i) / static_cast<double>(pts.size()));
            if (i > 0) CHECK(t[i] > t[i - 1]);
        }
    }

    // ---------- 9) 单参数曲线拟合：点在直线上时，拟合结果应为同一直线 ----------
    {
        // 直线 x = 1 + 0.5t, y = 2 - 1.5t（即 y = 5 - 3x），取 5 个等间隔 t 的点
        std::vector<curve::Point> pts;
        for (int i = 0; i < 5; ++i) {
            const double t = 0.5 * i;
            pts.push_back({1.0 + t, 2.0 - 3.0 * t});
        }
        const auto t = curve::parameterizeUniform(pts);
        const auto curvePts = curve::fitParametricCurve(pts, t, 1, 100);  // 1 次 = 直线
        CHECK(curvePts.size() == 100);
        for (const auto& p : curvePts)
            CHECK(std::abs(p.y - (5.0 - 3.0 * p.x)) < 1e-6);
    }

    // ---------- 10) 参数仿射变换不改变曲线几何：t 换成 t' = 2t+1 结果一致 ----------
    {
        const std::vector<curve::Point> pts = {
            {0.0, 0.0}, {1.0, 1.0}, {0.5, 2.0}, {-1.0, 1.0}};
        const auto t1 = curve::parameterizeUniform(pts);
        std::vector<double> t2(pts.size());
        for (std::size_t i = 0; i < t2.size(); ++i) t2[i] = 2.0 * static_cast<double>(i) + 1.0;

        const auto c1 = curve::fitParametricCurve(pts, t1, 2, 100);
        const auto c2 = curve::fitParametricCurve(pts, t2, 2, 100);
        CHECK(c1.size() == c2.size());
        // 仿射重参数化只是改变 t 的刻度：同一参数值对应同一二维点，逐点应相等
        for (std::size_t k = 0; k < c1.size(); ++k) {
            CHECK(std::abs(c1[k].x - c2[k].x) < 1e-6);
            CHECK(std::abs(c1[k].y - c2[k].y) < 1e-6);
        }
    }

    // ---------- 11) Foley 参数化：无符号下溢越界回归测试 ----------
    // 早期实现用 std::max(i - 2, 0u) 取 p0 下标，i=1 时 i-2 下溢成巨大数，
    // 直接 pts[巨大数] 越界崩溃。修复后应正常返回、严格递增且已归一化。
    {
        const std::vector<curve::Point> pts = {
            {0.0, 0.0}, {1.0, 2.0}, {2.0, 0.5}, {3.0, 1.8}, {4.0, 1.0}};
        const auto t = curve::parameterizeFoley(pts);
        CHECK(t.size() == pts.size());
        CHECK(t[0] == 0.0);
        for (std::size_t i = 1; i < t.size(); ++i) {
            CHECK(t[i] > t[i - 1]);        // 严格递增
            CHECK(t[i] == t[i]);           // 非 NaN
            CHECK(t[i] <= 1.0 + 1e-12);    // 已归一化到 [0,1]
        }

        // 最小点数（2 个）也要安全
        const std::vector<curve::Point> two = {{0.0, 0.0}, {1.0, 1.0}};
        const auto t2 = curve::parameterizeFoley(two);
        CHECK(t2.size() == 2);
        CHECK(t2[1] > t2[0]);

        // 全部重合的点：应返回全 0（不允许出现 NaN 或崩溃）
        const std::vector<curve::Point> same = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
        const auto ts = curve::parameterizeFoley(same);
        CHECK(ts.size() == 3);
        for (const double v : ts) CHECK(v == v && v == 0.0);

        // 端到端：Foley 参数化结果可直接用于单参数曲线拟合
        const auto c = curve::fitParametricCurve(pts, t, 2, 100);
        CHECK(c.size() == 100);
    }

    // ---------- 12) 1D k-means 撒点：单中心=均值，多中心升序且在范围内 ----------
    {
        const std::vector<double> xs = {0.0, 0.1, 0.2, 0.3};
        const auto c1 = curve::kMeansCenters1D(xs, 1);
        CHECK(c1.size() == 1);
        CHECK(std::abs(c1[0] - 0.15) < 1e-12);  // 均值

        const std::vector<double> xs2 = {0.0, 0.1, 0.9, 1.0};
        const auto c2 = curve::kMeansCenters1D(xs2, 2);
        CHECK(c2.size() == 2);
        CHECK(c2[0] < c2[1]);                 // 升序
        CHECK(c2[0] >= 0.0 && c2[1] <= 1.0);  // 在数据范围内

        // k > n 时自动限制到 n
        const auto c3 = curve::kMeansCenters1D(xs2, 100);
        CHECK(c3.size() == 4);
    }

    // ---------- 13) Gauss 最小二乘拟合：单中心能精确恢复单个高斯峰 ----------
    {
        const double c = 0.4, sigma = 0.2, A = 3.0;
        std::vector<curve::Point> pts;
        for (int i = 0; i <= 20; ++i) {
            const double x = 0.04 * i;  // 0 ~ 0.8
            pts.push_back({x, A * curve::gaussNodeValue(x, c, sigma)});
        }
        std::vector<double> centers;
        const auto coefs = curve::fitGaussLeastSquares(pts, sigma, 1, centers);
        CHECK(coefs.size() == 1);
        CHECK(centers.size() == 1);
        CHECK(std::abs(centers[0] - c) < 1e-6);  // k-means 单中心 = 均值 = c
        CHECK(std::abs(curve::evaluateGaussLeastSquares(coefs, centers, sigma, c) - A) < 1e-6);

        const auto fit = curve::gaussLeastSquares(pts, sigma, 1, 200);
        CHECK(fit.size() == 200);

        // 中心数超过点数时自动限制，不崩溃且返回有效采样
        const auto fit2 = curve::gaussLeastSquares(pts, sigma, 999, 100);
        CHECK(fit2.size() == 100);
    }

    // ---------- 14) Gauss 基参数曲线拟合：样本数正确、无 NaN、退化输入安全 ----------
    {
        const std::vector<curve::Point> pts = {
            {0.0, 0.0}, {1.0, 2.0}, {2.0, 0.5}, {3.0, 1.8}, {4.0, 1.0}};
        const auto t = curve::parameterizeUniform(pts);
        const auto c = curve::fitParametricCurveGauss(pts, t, 0.3, 5, 100);
        CHECK(c.size() == 100);
        for (const auto& p : c) {
            CHECK(p.x == p.x);  // 非 NaN
            CHECK(p.y == p.y);
        }

        // t 数量不匹配 / 中心数非法 → 返回空
        std::vector<double> badT(pts.size() + 1, 0.0);
        CHECK(curve::fitParametricCurveGauss(pts, badT, 0.3, 5, 100).empty());
        CHECK(curve::fitParametricCurveGauss(pts, t, 0.3, 0, 100).empty());
    }

    // ---------- 15) 参数型三次样条（弦长参数化）：过全部点、无 NaN ----------
    {
        // x 非单调的有序点列：标量三次样条处理不了，参数型可以
        const std::vector<curve::Point> pts = {
            {0.1, 0.2}, {0.6, 0.9}, {0.4, 0.3}, {0.9, 0.7}, {0.3, 0.8}};
        const int per = 20;
        const auto c = curve::fitParametricSpline(pts, curve::CubicSplineType::Natural, per);
        CHECK(c.size() == (pts.size() - 1) * static_cast<std::size_t>(per));
        for (const auto& p : c) {
            CHECK(p.x == p.x);  // 非 NaN
            CHECK(p.y == p.y);
        }
        // 段起点边界处的采样点应精确等于数据点（x(tᵢ)=xᵢ, y(tᵢ)=yᵢ）
        for (std::size_t i = 0; i < pts.size() - 1; ++i) {
            const auto& q = c[i * static_cast<std::size_t>(per)];
            CHECK(std::abs(q.x - pts[i].x) < 1e-6);
            CHECK(std::abs(q.y - pts[i].y) < 1e-6);
        }
        // 终点
        CHECK(std::abs(c.back().x - pts.back().x) < 1e-6);
        CHECK(std::abs(c.back().y - pts.back().y) < 1e-6);

        // 夹持边界同样过点
        const auto c2 = curve::fitParametricSpline(pts, curve::CubicSplineType::Clamped, per);
        CHECK(c2.size() == c.size());
        CHECK(std::abs(c2.back().x - pts.back().x) < 1e-6);
        CHECK(std::abs(c2.back().y - pts.back().y) < 1e-6);
    }

    // ---------- 16) 分段 Bezier（Catmull-Rom 控制点构造）----------
    {
        // n == 2：不需要创建虚拟点，B1/B2 直接取两端型值点 → 直线段
        const std::vector<curve::Point> two = {{0.2, 0.3}, {0.8, 0.9}};
        const auto segs2 = curve::buildBezierSegments(two);
        CHECK(segs2.size() == 1);
        CHECK(std::abs(segs2[0].control[0].x - two[0].x) < 1e-12);
        CHECK(std::abs(segs2[0].control[0].y - two[0].y) < 1e-12);
        CHECK(std::abs(segs2[0].control[1].x - two[1].x) < 1e-12);
        CHECK(std::abs(segs2[0].control[1].y - two[1].y) < 1e-12);
        // u=0.5 处应正好落在两端型值点的中点
        const auto mid = curve::evaluateBezierSegment(two, 0, segs2[0], 0.5);
        CHECK(std::abs(mid.x - 0.5) < 1e-12);
        CHECK(std::abs(mid.y - 0.6) < 1e-12);

        // 插值性质：每段 u=0 → P_i、u=1 → P_{i+1}，曲线过全部型值点
        const std::vector<curve::Point> pts = {
            {0.1, 0.2}, {0.6, 0.9}, {0.4, 0.3}, {0.9, 0.7}};
        const auto segs = curve::buildBezierSegments(pts);
        CHECK(segs.size() == pts.size() - 1);
        for (std::size_t i = 0; i < segs.size(); ++i) {
            const auto a = curve::evaluateBezierSegment(pts, static_cast<int>(i), segs[i], 0.0);
            const auto b = curve::evaluateBezierSegment(pts, static_cast<int>(i), segs[i], 1.0);
            CHECK(std::abs(a.x - pts[i].x) < 1e-12 && std::abs(a.y - pts[i].y) < 1e-12);
            CHECK(std::abs(b.x - pts[i + 1].x) < 1e-12 && std::abs(b.y - pts[i + 1].y) < 1e-12);
        }

        // 采样入口：样本数 = (n−1) × perSegment，且无 NaN
        const auto c = curve::bezierCatmullRomInterpolate(pts, 20);
        CHECK(c.size() == (pts.size() - 1) * 20);
        for (const auto& p : c) CHECK(p.x == p.x && p.y == p.y);

        // 均匀 Catmull-Rom 能精确还原二次抛物线 y = x² 的中间段（P1→P2）：
        // 该段 B1 = (4/3, 5/3)、B2 = (5/3, 8/3)，Bezier 求值应逐点落在抛物线上。
        //（边界段受镜像虚拟点影响不要求精确还原，这里只验中间段。）
        const std::vector<curve::Point> para = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0}};
        const auto segsP = curve::buildBezierSegments(para);
        CHECK(std::abs(segsP[1].control[0].x - 4.0 / 3.0) < 1e-12);
        CHECK(std::abs(segsP[1].control[0].y - 5.0 / 3.0) < 1e-12);
        CHECK(std::abs(segsP[1].control[1].x - 5.0 / 3.0) < 1e-12);
        CHECK(std::abs(segsP[1].control[1].y - 8.0 / 3.0) < 1e-12);
        for (int k = 0; k <= 20; ++k) {
            const double u = static_cast<double>(k) / 20.0;
            const auto p = curve::evaluateBezierSegment(para, 1, segsP[1], u);
            CHECK(std::abs(p.y - p.x * p.x) < 1e-9);
        }

        // 共线点：所有切线共线，曲线应严格保持在这条直线上
        const std::vector<curve::Point> line = {{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 1.0}};
        const auto lc = curve::bezierCatmullRomInterpolate(line, 25);
        for (const auto& p : lc) CHECK(std::abs(p.y - 1.0) < 1e-12);

        // 退化输入：少于 2 个点返回空
        CHECK(curve::bezierCatmullRomInterpolate({}).empty());
        CHECK(curve::bezierCatmullRomInterpolate({{0.0, 0.0}}).empty());
    }

    // ---------- 17) bezierCatmullRomInterpolate 支持外部传入段控制点 ----------
    {
        const std::vector<curve::Point> pts = {{0.0, 0.0}, {1.0, 2.0}, {2.0, 0.5}};
        // 手写"直线段"控制点：B1 = P_i、B2 = P_{i+1} → 每段几何上退化为线段
        std::vector<curve::BezierSegment> segs(2);
        segs[0].control[0] = pts[0];
        segs[0].control[1] = pts[1];
        segs[1].control[0] = pts[1];
        segs[1].control[1] = pts[2];

        const auto c = curve::bezierCatmullRomInterpolate(pts, segs, 20);
        CHECK(c.size() == 2 * 20);
        // 段 0（P0→P1）：几何上是 y = 2x 的线段
        for (int k = 0; k < 20; ++k) {
            const auto& p = c[static_cast<std::size_t>(k)];
            CHECK(std::abs(p.y - 2.0 * p.x) < 1e-9);
        }
        // 段 1（P1→P2）：几何上是 y = 3.5 − 1.5x 的线段
        for (int k = 0; k < 20; ++k) {
            const auto& p = c[static_cast<std::size_t>(20 + k)];
            CHECK(std::abs(p.y - (3.5 - 1.5 * p.x)) < 1e-9);
        }

        // 传入 buildBezierSegments 的结果应与不带段参数的版本完全一致
        const auto segsAuto = curve::buildBezierSegments(pts);
        const auto cAuto = curve::bezierCatmullRomInterpolate(pts, 20);
        const auto cGiven = curve::bezierCatmullRomInterpolate(pts, segsAuto, 20);
        CHECK(cAuto.size() == cGiven.size());
        for (std::size_t k = 0; k < cAuto.size(); ++k) {
            CHECK(std::abs(cAuto[k].x - cGiven[k].x) < 1e-12);
            CHECK(std::abs(cAuto[k].y - cGiven[k].y) < 1e-12);
        }

        // segs 数量不匹配（如 pts 变了但控制点没重建）→ 自动回退重建，仍有效
        CHECK(curve::bezierCatmullRomInterpolate(pts, {}, 20).size() == 40);
        std::vector<curve::BezierSegment> bad(3);  // n=3 需要 2 段
        CHECK(curve::bezierCatmullRomInterpolate(pts, bad, 20).size() == 40);
    }

    if (g_failures == 0) {
        std::printf("ALL TESTS PASSED\n");
        return 0;
    }
    std::printf("%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
