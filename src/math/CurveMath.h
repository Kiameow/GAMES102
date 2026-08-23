#pragma once

#include <Eigen/Dense>

#include <vector>

namespace curve
{

    // 二维数据点（纯数值结构，与 Qt 无关）
    struct Point
    {
        double x = 0.0;
        double y = 0.0;
    };

    double L2(const Point& p1, const Point& p2);
    double angleInRadians(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2);
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
    // Gauss 基函数最小二乘拟合（逼近，非插值）
    // 与插值版的区别：中心个数 m 可小于数据点数 n，中心位置由 1D k-means 对
    // x 坐标聚类得到（自适应数据分布），解过定最小二乘 min ||Φc − y|| 得到
    // 平滑逼近；m 和 σ 需搭配（σ 别远小于中心间距，否则设计矩阵列近共线、数值病态）。
    // ---------------------------------------------------------------------------

    // 1D k-means：把一维坐标 xs 聚类成 k 个中心，返回升序中心（k 自动限制到 [1,n]）。
    std::vector<double> kMeansCenters1D(const std::vector<double>& xs, int k);

    // 拟合：用 k-means 撒 numCenters 个高斯中心，解最小二乘，返回系数 c（长度 = 中心数）。
    // centers 输出实际使用的中心位置（升序）。sigma 为固定超参数（宽度）。
    std::vector<double> fitGaussLeastSquares(const std::vector<Point>& pts, double sigma,
        int numCenters, std::vector<double>& centers);

    double evaluateGaussLeastSquares(const std::vector<double>& coefs,
        const std::vector<double>& centers,
        double sigma, double x);

    // 一体式：k-means 撒点 + 拟合 + 采样，返回 samples 个点（直接给 PlotWidget 绘制）。
    std::vector<Point> gaussLeastSquares(const std::vector<Point>& pts, double sigma,
        int numCenters, int samples = 240);

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

    // ---------------------------------------------------------------------------
    // 单参数曲线拟合（作业3）
    // 思路：任意有序点列 P₀..Pₙ₋₁ 本身不是单值函数，但可以引入参数 t：
    //      把点列看成 (tᵢ, Pᵢ)，拆成两个标量点列 (tᵢ, xᵢ) 和 (tᵢ, yᵢ)，
    //      分别拟合 x(t)、y(t)，再合成二维曲线 (x(t), y(t))。
    // ---------------------------------------------------------------------------

    // 参数化：给定有序点列，计算每个点的参数值 tᵢ（约定 t₀ = 0，严格递增）。
    // 均匀参数化：tᵢ = i/n —— 与点之间的距离无关（最简单；点距不均匀时
    //             曲线在疏密变化处"速度"不均匀）。
    // 弦长参数化：tᵢ = tᵢ₋₁ + |Pᵢ − Pᵢ₋₁|，再归一化到 [0,1] —— 参数与几何弧长
    //             成正比，点密处 t 走得慢、点疏处走得快，曲线速度更均匀。
    // 中心参数化：tᵢ = tᵢ₋₁ + sqrt(|Pᵢ − Pᵢ₋₁|)，再归一化 —— 对急弯/尖角更友好，
    //             可缓解弦长法在曲率突变处容易出现的过冲/扭曲。
    // Foley 参数化：在弦长法基础上，按相邻两边的夹角给每段弦长加修正项
    //             （急弯处参数间隔被拉大），进一步抑制曲率突变处的扭曲。
    // 说明：曲线几何形状对参数做仿射变换 (t → a·t + b) 不变，因此是否归一化
    //      不影响最终曲线形状，归一化只是为了各方法数值尺度统一、便于对比。
    std::vector<double> parameterizeUniform(const std::vector<Point>& pts);
    std::vector<double> parameterizeChordal(const std::vector<Point>& pts);
    std::vector<double> parameterizeCentripetal(const std::vector<Point>& pts);
    std::vector<double> parameterizeFoley(const std::vector<Point>& pts);

    // 单参数曲线拟合：用参数化得到的 t 把点列拆成 (t,x)、(t,y) 两个标量点列，
    // 各自做 degree 次最小二乘多项式拟合，再在统一的 t 网格上采样合成二维曲线。
    // t 来自任意一种参数化方法（后续弦长/中心/Foley 参数化都直接喂进来即可）；
    // 返回 samples 个二维曲线点，可直接交给 PlotWidget 绘制。
    std::vector<Point> fitParametricCurve(const std::vector<Point>& pts,
        const std::vector<double>& t,
        int degree,
        int samples = 240);

    // 同 fitParametricCurve，但两个标量拟合改用 Gauss 基最小二乘
    //（k-means 撒 numCenters 个中心、宽度 sigma），供"参数拟合基函数"选 Gauss 时使用。
    std::vector<Point> fitParametricCurveGauss(const std::vector<Point>& pts,
        const std::vector<double>& t,
        double sigma,
        int numCenters,
        int samples = 240);

    // 三次样条插值（三转角法 / 斜率法）
    // fitParametricSpline 是通用入口：对任意有序点列先用弦长参数化得到 tᵢ
    //（严格递增），再对 (tᵢ,xᵢ)、(tᵢ,yᵢ) 分别做标量三次样条并合成二维曲线，
    // 因此 x、y 无需单调。连续性由标量样条保证：关于 t 的分量 x(t)、y(t)
    // 都是 C²，合成曲线也是 C²，参数化只影响形状。
    // Natural=自然边界(M₀=Mₙ=0)，Clamped=夹持边界（端点斜率用相邻差分估算）。
    // controls（可选，见 VertexTangent）：per-node 切线覆盖，被编辑节点不再
    // 要求 C²——Smooth 保持 C¹（共享切线）、Straight 保持 G¹（共线）、
    // Corner 允许 G⁰（拐角）。controls 为空 = 全部 Free，与旧行为一致。
    enum CubicSplineType
    {
        Natural,
        Clamped
    };

    // 顶点切线控制（样条交互式 G¹/G⁰ 编辑用，per-node 切线覆盖）：
    //  - Free:     不编辑，节点保持 C²（切线由系统解出）
    //  - Smooth:   t 为共享切线 dP/dt（左右一致 → 该点 C¹，C² 在该点断开）
    //  - Straight: tLeft/tRight 共线但长度可不同（G¹，共线由 UI 侧保证）
    //  - Corner:   tLeft/tRight 完全独立（G⁰，允许拐角）
    struct VertexTangent {
        enum class Mode { Free, Smooth, Straight, Corner };
        Mode mode = Mode::Free;
        Point t;        // Smooth：共享切线向量
        Point tLeft;    // Straight/Corner：左侧切线向量
        Point tRight;   // Straight/Corner：右侧切线向量
    };

    std::vector<Point> fitParametricSpline(const std::vector<Point>& pts,
                                           CubicSplineType splineType,
                                           int perSegmentSamples = 80,
                                           const std::vector<VertexTangent>& controls = {});

    // Hermite 分段插值求值：第 j 段 [x_j, x_{j+1}]（长 h），左端 (f_j, m_j^R)、
    // 右端 (f_{j+1}, m_{j+1}^L)，在 x 处取值（u=(x−x_j)/h ∈ [0,1]）。
    // 这是三转角法三次样条的段表示基础：整条样条 = 各段 Hermite 拼接。
    double hermiteEvaluate(double x, double xj, double h, double fj, double mjR,
                           double fj1, double mj1L);
    std::vector<double> cubicSplineCoefs(const std::vector<Point>& pts, CubicSplineType splineType);
    double evaluateCubicSplineSegment(const Point& stP, const std::vector<double>& coefs, double x);


    // Bezier Curve
    double binomialCoef(unsigned n, unsigned i);

    // 三次 Bernstein 基求值：pts 恰好 4 个控制点（B0..B3），参数 t ∈ [0,1]，
    // B(t) = Σ Bᵢ·C(3,i)·t^i·(1−t)^(3−i)。
    Point evaluateCubicBernsteinBasisCurve(const std::vector<Point>& pts, double t);

    // ---------------------------------------------------------------------------
    // 分段三次 Bezier（Catmull-Rom 控制点构造）
    // 每段 [P_i, P_{i+1}] 是一条三次 Bezier，4 个控制点为：
    //     B0 = P_i                          （型值点 = 段起点）
    //     B1 = P_i + (P_{i+1} − P_{i−1})/6  （由 Catmull-Rom 切线推出的中间控制点）
    //     B2 = P_{i+1} − (P_{i+2} − P_i)/6  （由 Catmull-Rom 切线推出的中间控制点）
    //     B3 = P_{i+1}                      （型值点 = 段终点）
    // 边界段需要虚拟点 P_{−1}、P_n（catmullRomVirtualPoint 端点镜像外推）；
    // n == 2 时只有一段、没有相邻点，不需要创建虚拟点，B1/B2 直接取 P0/P1（直线段）。
    // 曲线过全部型值点，且在型值点处切线连续（C¹，由 Catmull-Rom 共享切线保证）。
    // ---------------------------------------------------------------------------

    // 每段的统一控制点结构：管理本段的两个中间控制点（B1、B2）。
    // 段端点 B0、B3 就是两侧型值点（pts[i]、pts[i+1]），由调用方持有，不重复存储。
    struct BezierSegment {
        Point control[2];  // control[0] = B1（起点侧）、control[1] = B2（终点侧）
    };

    // Catmull-Rom 虚拟点：下标越界（-1 / n）时用端点镜像外推生成虚拟点
    // （P_{−1} = 2·P₀ − P₁，P_n = 2·P_{n−1} − P_{n−2}），使边界段也能套用统一切线公式；
    // 下标在范围内时原样返回数据点。
    Point catmullRomVirtualPoint(const std::vector<Point>& pts, int index);

    // 构造第 i 段（P_i → P_{i+1}）的中间控制点；n == 2 时直接返回 B1=P0、B2=P1。
    BezierSegment buildBezierSegment(const std::vector<Point>& pts, int i);

    // 一次生成全部 n−1 段的中间控制点（每段一个 BezierSegment）。
    std::vector<BezierSegment> buildBezierSegments(const std::vector<Point>& pts);

    // 第 i 段在参数 u ∈ [0,1] 处求值：结合本段的 control[] 与自身两个型值点
    // pts[i]、pts[i+1] 做三次 Bernstein 插值（u 越界自动 clamp 到 [0,1]）。
    Point evaluateBezierSegment(const std::vector<Point>& pts, int i,
                                const BezierSegment& seg, double u);

    // 一体式入口：Catmull-Rom 控制点 + 每段 perSegmentSamples 个采样点，
    // 返回可直接交给 PlotWidget 绘制的点列。
    std::vector<Point> bezierCatmullRomInterpolate(const std::vector<Point>& pts,
                                                   int perSegmentSamples = 80);

    // 同上一体式入口，但使用调用方给定的段控制点 segs（支持 GUI 手动拖过控制点：
    // 拖过的 B1/B2 直接生效，而不是从 pts 重新按 Catmull-Rom 推导）。
    // segs 长度必须等于 n−1，否则自动回退到 buildBezierSegments(pts) 重建。
    std::vector<Point> bezierCatmullRomInterpolate(const std::vector<Point>& pts,
                                                   const std::vector<BezierSegment>& segs,
                                                   int perSegmentSamples = 80);
}  // namespace curve
