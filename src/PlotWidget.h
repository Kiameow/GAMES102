#pragma once

#include "math/CurveMath.h"

#include <QColor>
#include <QPoint>
#include <QPointF>
#include <QPolygonF>
#include <QRectF>
#include <QString>
#include <QStringList>
#include <QVector>
#include <QWidget>

#include <functional>
#include <vector>

class QTimer;

// 一个算法层的"计算函数"：输入数据点、输出采样点（与 CurveMath 的算法签名一致）
using CurveCompute =
    std::function<std::vector<curve::Point>(const std::vector<curve::Point>&, int samples)>;

// 一条曲线算法层：算法 = 名字 + 颜色 + 计算函数。
// 注册后由 PlotWidget 统一管理：开关可见、重算采样、按颜色绘制、图例。
struct CurveLayer {
    QString key;    // 唯一标识（按钮 / 状态栏用它定位）
    QString name;   // 显示名（按钮、图例）
    QColor color;   // 曲线颜色
    bool visible = false;
    CurveCompute compute;   // 数学算法入口（外部结果层可以为空）
    QPolygonF samples;      // 最近一次计算的采样点（数据坐标系，paintEvent 里绘制）
    bool external = false;  // 结果来自外部（如 Python 训练），rebuildCurve 不重算它
};

// 参数曲线拟合（作业3）的基函数选择：幂基多项式 / Gauss 基（k-means 撒点）
enum class ParametricBasis { Power, Gauss };

// 绘图画布：
//   - 数据空间固定为 [0,1] x [0,1]（左下角为原点）
//   - 鼠标左键添加红点，右键删除最近的红点
//   - 任意多条算法曲线可同时显示（"显示所有"），各自颜色 + 图例
class PlotWidget : public QWidget {
    Q_OBJECT
public:
    explicit PlotWidget(QWidget* parent = nullptr);

    // ---- 算法层注册（通用入口） ----
    void registerLayer(const CurveLayer& layer);  // key 相同则更新，保留可见状态
    void setLayerVisible(const QString& key, bool visible);
    bool isLayerVisible(const QString& key) const;
    void setAllLayersVisible(bool visible);
    QStringList layerKeys() const;
    QString layerName(const QString& key) const;

    // ---- 数据点 ----
    int pointCount() const { return m_points.size(); }
    void clearPoints();
    const QVector<QPointF>& points() const { return m_points; }

    // ---- 外部结果回填（如 Python 训练结果；数据点变化时自动作废） ----
    void setExternalCurve(const QString& key, std::vector<curve::Point> samples);

    // ---- 最小二乘次数（逼近类算法通过 approximationDegree() 读取） ----
    int approximationDegree() const { return m_degree; }
    void setApproximationDegree(int degree);

    // ---- 高斯基函数宽度 σ（Gauss 类算法通过 gaussianSigma() 读取） ----
    double gaussianSigma() const { return m_sigma; }
    void setGaussianSigma(double sigma);

    // ---- Gauss 基拟合中心个数（Gauss 最小二乘算法通过 gaussCenters() 读取） ----
    int gaussCenters() const { return m_gaussCenters; }
    void setGaussCenters(int centers);

    // ---- 参数曲线拟合的基函数（4 个"参数曲线-*"按钮共用） ----
    ParametricBasis parametricBasis() const { return m_paramBasis; }
    void setParametricBasis(ParametricBasis basis);

    // ---- 岭回归正则系数 λ（Ridge 类算法通过 lambda() 读取） ----
    double lambda() const { return m_lambda; }
    void setLambda(double lambda);

    // ---- 顶点切线编辑（右键选中节点，拖切线柄端点调整连续性） ----
    int selectedIndex() const { return m_selectedIndex; }
    curve::VertexTangent::Mode selectedVertexMode() const;  // 选中节点的模式（无选中返回 pending）
    void setSelectedVertexMode(curve::VertexTangent::Mode mode);  // 应用到选中节点并记住 pending
    const std::vector<curve::VertexTangent>& tangentControls() const { return m_tangentControls; }
    void removeLastPoint();  // 删除最近放置的点（替代原右键删除）

    // ---- 分段 Bezier 的中间控制点（全局管理，类似 m_points；可左键拖动） ----
    // 每段两个（B1、B2），与 m_points 平行：pts 变化时自动按 Catmull-Rom 重建，
    // 但用户手动拖过的控制点（pinned）会保留原位。曲线计算入口读取该列表。
    const std::vector<curve::BezierSegment>& bezierSegments() const { return m_bezierSegments; }

signals:
    void pointsChanged(int count);
    void nodeSelected(int index);  // 右键选中节点（-1 = 取消选中）

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;

private slots:
    void onDragTimeout();  // 长按 0.5s 计时到：进入点拖动模式

private:
    QRectF dataRect() const;
    QPointF dataToWidget(const QPointF& data) const;
    QPointF widgetToData(const QPointF& widget) const;
    QPolygonF toWidgetCoords(const QPolygonF& data) const;
    int nearestPointIndex(const QPointF& widgetPos, double maxDistPx) const;
    void rebuildCurve();
    void clearExternalResults();  // 数据点变化时作废全部外部结果
    void cancelDragArm();         // 取消"等待长按"状态（快速点击或移开）
    void selectNode(int index);   // 右键选中节点（套用 pending 模式到新节点）
    void applyVertexModeTo(int index, curve::VertexTangent::Mode mode);
    curve::Point estimateTangent(int index) const;  // 用中央差分估计节点切线 dP/dt
    bool handleEndpoints(QPointF& rightEnd, QPointF& leftEnd) const;  // 选中节点柄端点(widget坐标)
    QPointF dataVecToWidget(const QPointF& v) const;  // 数据向量 → widget 向量（y 翻转）
    void straightenLeft(curve::VertexTangent& c);   // 直线顶点：右侧拖动后左侧共线
    void straightenRight(curve::VertexTangent& c);  // 直线顶点：左侧拖动后右侧共线

    // ---- 分段 Bezier 控制点辅助 ----
    void rebuildBezierControls();  // pts 变化时重建控制点（保留用户拖过的 pinned 覆盖）
    int nearestBezierIndex(const QPointF& widgetPos, double maxDistPx) const;  // 命中控制点(全局下标,-1=无)
    QPointF bezierControlWidgetPos(int flat) const;  // 控制点全局下标(2i+k) → widget 坐标

    QVector<QPointF> m_points;
    std::vector<CurveLayer> m_layers;  // 已注册的算法层
    int m_degree = 3;
    double m_sigma = 0.3;
    double m_lambda = 0.1;
    int m_gaussCenters = 5;
    ParametricBasis m_paramBasis = ParametricBasis::Power;

    // ---- 点拖动状态（左键长按 0.5s 选中并拖动型值点） ----
    QTimer* m_dragTimer = nullptr;  // 长按计时器（单次 500ms）
    int m_armIndex = -1;            // 按下时命中的点，等待 0.5s 进入拖动
    int m_dragIndex = -1;           // 正在拖动的点（-1 = 无）
    QPoint m_pressWidgetPos;        // 按下时的位置（用于移开取消误触）

    // ---- 顶点切线编辑状态 ----
    int m_selectedIndex = -1;  // 右键选中的节点（-1 = 无）
    std::vector<curve::VertexTangent> m_tangentControls;  // 与 m_points 平行的切线控制
    int m_dragHandle = 0;       // 切线柄拖动：+1 右柄 / -1 左柄 / 0 无
    curve::VertexTangent::Mode m_pendingVertexMode =
        curve::VertexTangent::Mode::Smooth;  // 右键选中新节点时套用的模式

    // ---- 分段 Bezier 控制点状态（与 m_points 平行，全局管理） ----
    std::vector<curve::BezierSegment> m_bezierSegments;  // 当前每段的中间控制点（含用户拖动覆盖）
    std::vector<bool> m_bezierPinned;  // 与 2×(n−1) 个控制点平行：用户手动拖过 = true（pts 变化时保留）
    int m_dragBezier = -1;             // 正在拖动的控制点全局下标（2×段 + k；-1 = 无）

    static constexpr double kMargin = 30.0;  // 四周留白（像素）
};
