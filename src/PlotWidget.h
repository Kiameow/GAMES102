#pragma once

#include "math/CurveMath.h"

#include <QColor>
#include <QPointF>
#include <QPolygonF>
#include <QRectF>
#include <QString>
#include <QStringList>
#include <QVector>
#include <QWidget>

#include <functional>
#include <vector>

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

    // ---- 岭回归正则系数 λ（Ridge 类算法通过 lambda() 读取） ----
    double lambda() const { return m_lambda; }
    void setLambda(double lambda);

signals:
    void pointsChanged(int count);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;

private:
    QRectF dataRect() const;
    QPointF dataToWidget(const QPointF& data) const;
    QPointF widgetToData(const QPointF& widget) const;
    QPolygonF toWidgetCoords(const QPolygonF& data) const;
    int nearestPointIndex(const QPointF& widgetPos, double maxDistPx) const;
    void rebuildCurve();
    void clearExternalResults();  // 数据点变化时作废全部外部结果

    QVector<QPointF> m_points;
    std::vector<CurveLayer> m_layers;  // 已注册的算法层
    int m_degree = 3;
    double m_sigma = 0.3;
    double m_lambda = 0.1;

    static constexpr double kMargin = 30.0;  // 四周留白（像素）
};
