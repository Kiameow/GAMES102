#include "PlotWidget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>

#include <algorithm>

namespace {
constexpr double kMinAddDistancePx = 8.0;    // 添加点时避免与已有红点重叠
constexpr double kRemoveDistancePx = 20.0;   // 右键删除的判定半径
constexpr int kCurveSamples = 240;           // 曲线采样点数
}  // namespace

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent) {
    setCursor(Qt::CrossCursor);
    setMinimumSize(480, 360);
}

QRectF PlotWidget::dataRect() const {
    const double w = std::max(1.0, width() - 2.0 * kMargin);
    const double h = std::max(1.0, height() - 2.0 * kMargin);
    return QRectF(kMargin, kMargin, w, h);
}

QPointF PlotWidget::dataToWidget(const QPointF& d) const {
    const QRectF r = dataRect();
    return QPointF(r.left() + d.x() * r.width(),
                   r.bottom() - d.y() * r.height());  // y 轴翻转
}

QPointF PlotWidget::widgetToData(const QPointF& w) const {
    const QRectF r = dataRect();
    const double x = (w.x() - r.left()) / r.width();
    const double y = (r.bottom() - w.y()) / r.height();
    return QPointF(std::clamp(x, 0.0, 1.0), std::clamp(y, 0.0, 1.0));
}

QPolygonF PlotWidget::toWidgetCoords(const QPolygonF& data) const {
    QPolygonF out;
    out.reserve(data.size());
    for (const QPointF& p : data) out << dataToWidget(p);
    return out;
}

int PlotWidget::nearestPointIndex(const QPointF& widgetPos, double maxDistPx) const {
    int best = -1;
    double bestDist = maxDistPx * maxDistPx;
    for (int i = 0; i < m_points.size(); ++i) {
        const QPointF w = dataToWidget(m_points[i]);
        const double dx = w.x() - widgetPos.x();
        const double dy = w.y() - widgetPos.y();
        const double d2 = dx * dx + dy * dy;
        if (d2 <= bestDist) {
            bestDist = d2;
            best = i;
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// 算法层管理
// ---------------------------------------------------------------------------

void PlotWidget::registerLayer(const CurveLayer& layer) {
    for (CurveLayer& l : m_layers) {
        if (l.key == layer.key) {
            // 同名层更新（名字/颜色/计算函数），保留可见状态
            l.name = layer.name;
            l.color = layer.color;
            l.compute = layer.compute;
            if (l.visible) rebuildCurve();
            return;
        }
    }
    m_layers.push_back(layer);  // 新层默认不可见
}

void PlotWidget::setLayerVisible(const QString& key, bool visible) {
    for (CurveLayer& l : m_layers) {
        if (l.key == key) {
            if (l.visible == visible) return;
            l.visible = visible;
            rebuildCurve();
            return;
        }
    }
}

bool PlotWidget::isLayerVisible(const QString& key) const {
    for (const CurveLayer& l : m_layers)
        if (l.key == key) return l.visible;
    return false;
}

void PlotWidget::setAllLayersVisible(bool visible) {
    bool changed = false;
    for (CurveLayer& l : m_layers) {
        if (l.visible != visible) {
            l.visible = visible;
            changed = true;
        }
    }
    if (changed) rebuildCurve();
}

QStringList PlotWidget::layerKeys() const {
    QStringList keys;
    for (const CurveLayer& l : m_layers) keys << l.key;
    return keys;
}

QString PlotWidget::layerName(const QString& key) const {
    for (const CurveLayer& l : m_layers)
        if (l.key == key) return l.name;
    return QString();
}

// 把外部（如 Python）算好的采样点直接写入某层并显示
void PlotWidget::setExternalCurve(const QString& key, std::vector<curve::Point> samples) {
    for (CurveLayer& l : m_layers) {
        if (l.key == key) {
            l.samples.clear();
            for (const auto& s : samples) l.samples << QPointF(s.x, s.y);
            l.external = true;
            l.visible = true;  // 回填成功即显示
            update();
            return;
        }
    }
}

// ---------------------------------------------------------------------------
// 数据点 / 曲线重算
// ---------------------------------------------------------------------------

void PlotWidget::setApproximationDegree(int degree) {
    m_degree = std::max(1, degree);
    rebuildCurve();  // 次数可能影响可见的逼近曲线，一律重算
}

void PlotWidget::setGaussianSigma(double sigma) {
    m_sigma = sigma;
    rebuildCurve();  // σ 可能影响可见的高斯曲线，一律重算
}

void PlotWidget::setLambda(double lambda) {
    m_lambda = lambda;
    rebuildCurve();  // λ 可能影响可见的岭回归曲线，一律重算
}

void PlotWidget::clearPoints() {
    m_points.clear();
    clearExternalResults();
    rebuildCurve();
    emit pointsChanged(0);
}

void PlotWidget::clearExternalResults() {
    bool changed = false;
    for (CurveLayer& l : m_layers) {
        if (l.external) {
            l.samples.clear();
            l.external = false;
            changed = true;
        }
    }
    if (changed) update();
}

void PlotWidget::rebuildCurve() {
    std::vector<curve::Point> pts;
    pts.reserve(m_points.size());
    for (const QPointF& p : m_points) pts.push_back({p.x(), p.y()});

    // 逐个可见层调用算法，得到采样点；外部结果层不重算
    for (CurveLayer& l : m_layers) {
        if (l.external) continue;  // 结果由外部回填（如 Python），数据变化时由 clearExternalResults 作废
        l.samples.clear();
        if (!l.visible || pts.size() < 2 || !l.compute) continue;
        const auto samples = l.compute(pts, kCurveSamples);
        for (const auto& s : samples) l.samples << QPointF(s.x, s.y);
    }
    update();
}

void PlotWidget::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        // 与已有红点距离过近时忽略，避免堆叠
        if (nearestPointIndex(event->pos(), kMinAddDistancePx) < 0) {
            m_points.append(widgetToData(event->pos()));
            clearExternalResults();  // 数据变了，外部结果作废
            rebuildCurve();
            emit pointsChanged(m_points.size());
        }
    } else if (event->button() == Qt::RightButton) {
        const int idx = nearestPointIndex(event->pos(), kRemoveDistancePx);
        if (idx >= 0) {
            m_points.removeAt(idx);
            clearExternalResults();
            rebuildCurve();
            emit pointsChanged(m_points.size());
        }
    }
    QWidget::mousePressEvent(event);
}

void PlotWidget::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.fillRect(rect(), Qt::white);

    const QRectF dr = dataRect();

    // 0.1 间距的浅色网格
    painter.setPen(QPen(QColor(0xe8, 0xe8, 0xe8), 1));
    for (int i = 0; i <= 10; ++i) {
        const double fx = dr.left() + dr.width() * i / 10.0;
        painter.drawLine(QPointF(fx, dr.top()), QPointF(fx, dr.bottom()));
        const double fy = dr.top() + dr.height() * i / 10.0;
        painter.drawLine(QPointF(dr.left(), fy), QPointF(dr.right(), fy));
    }

    // 坐标框
    painter.setPen(QPen(Qt::black, 1.5));
    painter.drawRect(dr);

    // 各算法曲线（按各自颜色，可同时显示）
    for (const CurveLayer& l : m_layers) {
        if (!l.visible || l.samples.isEmpty()) continue;
        painter.setPen(QPen(l.color, 2));
        painter.drawPolyline(toWidgetCoords(l.samples));
    }

    // 数据点（红）
    painter.setPen(QPen(Qt::black, 1));
    painter.setBrush(Qt::red);
    for (const QPointF& p : m_points) {
        painter.drawEllipse(dataToWidget(p), 4.0, 4.0);
    }

    // 顶部提示 + 图例
    QFont f = font();
    f.setPointSize(10);
    painter.setFont(f);
    const QRectF hint(dr.left(), 6, dr.width(), 20);
    painter.setPen(QColor(0x88, 0x88, 0x88));
    painter.drawText(hint, Qt::AlignLeft | Qt::AlignVCenter,
                     QStringLiteral("左键添加红点 · 右键删除最近点 · 当前 %1 个点")
                         .arg(m_points.size()));

    // 图例：每个可见算法一行（圆角色块 + 名字），行距留宽避免拥挤
    int legendY = 30;
    const int rowH = 22;
    for (const CurveLayer& l : m_layers) {
        if (!l.visible) continue;
        painter.setPen(Qt::NoPen);
        painter.setBrush(l.color);
        painter.drawRoundedRect(QRectF(dr.left(), legendY + 3, 16, 13), 3, 3);
        painter.setPen(l.color);
        painter.drawText(QRectF(dr.left() + 24, legendY, dr.width() - 24, rowH),
                         Qt::AlignLeft | Qt::AlignVCenter, l.name);
        legendY += rowH;
    }
}
