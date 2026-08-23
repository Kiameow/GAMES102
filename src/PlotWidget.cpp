#include "PlotWidget.h"

#include <QGuiApplication>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QTimer>

#include <algorithm>
#include <cmath>

namespace {
constexpr double kMinAddDistancePx = 8.0;    // 添加点时避免与已有红点重叠
constexpr double kRemoveDistancePx = 20.0;   // 右键删除的判定半径
constexpr double kBezierHitRadiusPx = 10.0;  // Bezier 控制点（靛蓝方块）的命中半径
constexpr double kDragArmMoveCancelPx = 15.0; // 等待长按时移动超过该值则取消（防误触）
constexpr int kDragArmTimeoutMs = 500;        // 长按进入拖动的阈值
constexpr int kCurveSamples = 240;           // 曲线采样点数
// 切线柄显示缩放：存储值 = 真实一阶导 dP/dt（数学用它），显示值 = 存储值 × 该系数
//（减小可让柄更短；拖动时端点仍跟手，见 mouseMoveEvent 里的 1/系数 换算）
constexpr double kHandleDisplayScale = 0.2;
}  // namespace

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent) {
    setCursor(Qt::CrossCursor);
    setMinimumSize(480, 360);
    m_dragTimer = new QTimer(this);
    m_dragTimer->setSingleShot(true);
    m_dragTimer->setInterval(kDragArmTimeoutMs);
    connect(m_dragTimer, &QTimer::timeout, this, &PlotWidget::onDragTimeout);
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

void PlotWidget::setGaussCenters(int centers) {
    m_gaussCenters = std::max(1, centers);
    rebuildCurve();  // 中心数影响可见的 Gauss 拟合曲线
}

void PlotWidget::setParametricBasis(ParametricBasis basis) {
    if (m_paramBasis == basis) return;
    m_paramBasis = basis;
    rebuildCurve();  // 基函数变化影响可见的参数曲线
}

void PlotWidget::setLambda(double lambda) {
    m_lambda = lambda;
    rebuildCurve();  // λ 可能影响可见的岭回归曲线，一律重算
}

void PlotWidget::clearPoints() {
    m_points.clear();
    m_tangentControls.clear();
    m_selectedIndex = -1;
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
    // 先把分段 Bezier 的中间控制点与当前 m_points 对齐（用户拖过的控制点保留原位）
    rebuildBezierControls();

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

// ---------------------------------------------------------------------------
// 分段 Bezier 控制点（全局管理，类似 m_points；可左键拖动，无右键切线）
// ---------------------------------------------------------------------------

// pts 变化后重建全部段的中间控制点：先按 Catmull-Rom 从 m_points 推导，
// 再把"用户手动拖过"（pinned）的控制点原位覆盖回去（段索引仍有效时）。
void PlotWidget::rebuildBezierControls() {
    const int n = static_cast<int>(m_points.size());
    if (n < 2) {
        m_bezierSegments.clear();
        m_bezierPinned.clear();
        m_dragBezier = -1;
        return;
    }

    // 备份旧的 pinned 覆盖（位置 + 标记），重建后再按段/控制点下标恢复
    const auto oldSegs = m_bezierSegments;
    const auto oldPinned = m_bezierPinned;

    std::vector<curve::Point> pts;
    pts.reserve(static_cast<std::size_t>(n));
    for (const QPointF& p : m_points) pts.push_back({p.x(), p.y()});
    m_bezierSegments = curve::buildBezierSegments(pts);
    m_bezierPinned.assign(2 * static_cast<std::size_t>(n - 1), false);

    const int oldSegCount = static_cast<int>(oldSegs.size());
    for (int i = 0; i < n - 1 && i < oldSegCount; ++i) {
        for (int k = 0; k < 2; ++k) {
            const std::size_t flat = 2 * static_cast<std::size_t>(i) + static_cast<std::size_t>(k);
            if (flat < oldPinned.size() && oldPinned[flat]) {
                m_bezierSegments[static_cast<std::size_t>(i)].control[k] =
                    oldSegs[static_cast<std::size_t>(i)].control[k];
                m_bezierPinned[flat] = true;
            }
        }
    }

    // 点数减少导致拖动的控制点不存在了 → 取消拖动
    if (m_dragBezier >= 2 * (n - 1)) m_dragBezier = -1;
}

// 命中 Bezier 中间控制点：返回全局下标 2×段+k（-1 = 无）
int PlotWidget::nearestBezierIndex(const QPointF& widgetPos, double maxDistPx) const {
    int best = -1;
    double bestDist = maxDistPx * maxDistPx;
    for (std::size_t s = 0; s < m_bezierSegments.size(); ++s) {
        for (int k = 0; k < 2; ++k) {
            const QPointF w = dataToWidget(QPointF(m_bezierSegments[s].control[k].x,
                                                   m_bezierSegments[s].control[k].y));
            const double dx = w.x() - widgetPos.x();
            const double dy = w.y() - widgetPos.y();
            const double d2 = dx * dx + dy * dy;
            if (d2 <= bestDist) {
                bestDist = d2;
                best = static_cast<int>(2 * s + k);
            }
        }
    }
    return best;
}

// 控制点全局下标（2×段+k）→ widget 坐标（越界返回空点）
QPointF PlotWidget::bezierControlWidgetPos(int flat) const {
    const int s = flat / 2;
    const int k = flat % 2;
    if (s < 0 || s >= static_cast<int>(m_bezierSegments.size())) return QPointF();
    return dataToWidget(QPointF(m_bezierSegments[static_cast<std::size_t>(s)].control[k].x,
                                m_bezierSegments[static_cast<std::size_t>(s)].control[k].y));
}

void PlotWidget::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        // 0) 切线柄端点优先命中（选中节点有柄时），立即进入柄拖动
        QPointF rightEnd, leftEnd;
        if (handleEndpoints(rightEnd, leftEnd)) {
            const double dR = std::hypot(rightEnd.x() - event->pos().x(),
                                         rightEnd.y() - event->pos().y());
            const double dL = std::hypot(leftEnd.x() - event->pos().x(),
                                         leftEnd.y() - event->pos().y());
            if (dR <= kRemoveDistancePx || dL <= kRemoveDistancePx) {
                m_dragHandle = (dR <= dL) ? 1 : -1;
                update();
                return;
            }
        }
        // 0.5) Bezier 中间控制点：命中即立即进入拖动（无需长按；没有右键切线控制）。
        //      仅当该层勾选（方块可见）时参与命中；且只有比数据点更近时才优先
        //      （n==2 时 B1/P0、B2/P1 完全重合，数据点拖动优先，避免抢走型值点的长按拖点）。
        if (isLayerVisible(QStringLiteral("bezier_catmull_rom"))) {
            const int bz = nearestBezierIndex(event->pos(), kBezierHitRadiusPx);
            if (bz >= 0) {
                const int dpHit = nearestPointIndex(event->pos(), kMinAddDistancePx);
                bool preferBezier = dpHit < 0;
                if (dpHit >= 0) {
                    const QPointF bzW = bezierControlWidgetPos(bz);
                    const QPointF dpW = dataToWidget(m_points[dpHit]);
                    preferBezier = (bzW - event->pos()).manhattanLength() <
                                   (dpW - event->pos()).manhattanLength();
                }
                if (preferBezier) {
                    m_dragBezier = bz;
                    grabMouse(Qt::CrossCursor);
                    update();
                    return;
                }
            }
        }
        // 1) 命中已有红点：不立即加点，启动长按计时准备拖动
        const int hit = nearestPointIndex(event->pos(), kMinAddDistancePx);
        if (hit >= 0) {
            m_armIndex = hit;
            m_pressWidgetPos = event->pos();
            m_dragTimer->start();
            return;  // 等待 0.5s 或松开，不处理为"加点"
        }
        // 2) 空处：照旧添加红点（同步一条默认切线控制）
        m_points.append(widgetToData(event->pos()));
        m_tangentControls.push_back(curve::VertexTangent{});
        clearExternalResults();  // 数据变了，外部结果作废
        rebuildCurve();
        emit pointsChanged(m_points.size());
    } else if (event->button() == Qt::RightButton) {
        cancelDragArm();  // 右键期间取消等待，避免误入拖动
        const int idx = nearestPointIndex(event->pos(), kRemoveDistancePx);
        selectNode(idx);  // 命中 → 选中并出现切线柄；空处 → 取消选中
        return;
    }
    QWidget::mousePressEvent(event);
}

void PlotWidget::mouseMoveEvent(QMouseEvent* event) {
    // Bezier 控制点拖动：更新位置并标记 pinned（用户覆盖），实时重算曲线
    if (m_dragBezier >= 0) {
        const int s = m_dragBezier / 2;
        const int k = m_dragBezier % 2;
        if (s >= 0 && s < static_cast<int>(m_bezierSegments.size())) {
            const QPointF q = widgetToData(event->pos());
            m_bezierSegments[static_cast<std::size_t>(s)].control[k] = {q.x(), q.y()};
            m_bezierPinned[static_cast<std::size_t>(m_dragBezier)] = true;
            rebuildCurve();
        }
        return;
    }
    // 切线柄拖动：改选中节点的切线（实时重算曲线）
    if (m_dragHandle != 0 && m_selectedIndex >= 0 &&
        m_selectedIndex < static_cast<int>(m_points.size())) {
        const QPointF pkData = m_points[m_selectedIndex];
        const QPointF qData = widgetToData(event->pos());
        const QPointF vec = qData - pkData;  // 数据坐标下的切线向量（显示缩放的 1/s 换算）
        curve::VertexTangent& c = m_tangentControls[static_cast<std::size_t>(m_selectedIndex)];
        const double invScale = 1.0 / kHandleDisplayScale;
        if (c.mode == curve::VertexTangent::Mode::Smooth) {
            // 平滑顶点：拖任一端都设置共享切线（左端 = -vec）
            c.t = (m_dragHandle > 0)
                      ? curve::Point{vec.x() * invScale, vec.y() * invScale}
                      : curve::Point{-vec.x() * invScale, -vec.y() * invScale};
        } else if (c.mode == curve::VertexTangent::Mode::Straight ||
                   c.mode == curve::VertexTangent::Mode::Corner) {
            if (m_dragHandle > 0) {
                c.tRight = {vec.x() * invScale, vec.y() * invScale};
                if (c.mode == curve::VertexTangent::Mode::Straight) straightenLeft(c);
            } else {
                c.tLeft = {-vec.x() * invScale, -vec.y() * invScale};
                if (c.mode == curve::VertexTangent::Mode::Straight) straightenRight(c);
            }
        }
        rebuildCurve();
        return;
    }
    if (m_dragIndex >= 0) {
        // 拖动中：更新型值点位置，实时重算曲线（保持 C²）
        m_points[m_dragIndex] = widgetToData(event->pos());
        rebuildCurve();
        return;
    }
    // 等待长按期间大幅移动 → 取消（用户可能是误按）
    if (m_armIndex >= 0 &&
        (event->pos() - m_pressWidgetPos).manhattanLength() > kDragArmMoveCancelPx) {
        cancelDragArm();
    }
    QWidget::mouseMoveEvent(event);
}

void PlotWidget::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        if (m_dragBezier >= 0) {
            // Bezier 控制点拖动结束：位置已实时更新（pinned 保留），无需再重算
            m_dragBezier = -1;
            releaseMouse();
            rebuildCurve();
            return;
        }
        if (m_dragHandle != 0) {
            // 切线柄拖动结束：切线已实时更新，无需再重算
            m_dragHandle = 0;
            update();
            return;
        }
        if (m_dragIndex >= 0) {
            // 点拖动结束：固定该点位置
            m_dragIndex = -1;
            releaseMouse();
            clearExternalResults();  // 数据变了，外部结果作废
            rebuildCurve();
            emit pointsChanged(m_points.size());
            return;
        }
        // 快速点击（<0.5s）：取消等待，不产生任何改动
        cancelDragArm();
    }
    QWidget::mouseReleaseEvent(event);
}

void PlotWidget::onDragTimeout() {
    // 0.5s 到且左键仍按住、命中点仍存在 → 进入拖动模式
    if (!(QGuiApplication::mouseButtons() & Qt::LeftButton)) {
        m_armIndex = -1;
        return;
    }
    if (m_armIndex < 0 || m_armIndex >= m_points.size()) return;
    m_dragIndex = m_armIndex;
    m_armIndex = -1;
    grabMouse(Qt::CrossCursor);  // 拖动期间即使移出窗口也继续接收鼠标事件
    update();  // 绘制选中描边
}

void PlotWidget::cancelDragArm() {
    m_armIndex = -1;
    if (m_dragTimer) m_dragTimer->stop();
}

void PlotWidget::removeLastPoint() {
    if (m_points.isEmpty()) return;
    m_points.removeLast();
    if (!m_tangentControls.empty()) m_tangentControls.pop_back();
    if (m_selectedIndex >= static_cast<int>(m_points.size())) m_selectedIndex = -1;
    if (m_armIndex >= static_cast<int>(m_points.size())) m_armIndex = -1;
    if (m_dragIndex >= static_cast<int>(m_points.size())) m_dragIndex = -1;
    if (m_dragHandle != 0) {
        m_dragHandle = 0;
        releaseMouse();
    }
    clearExternalResults();
    rebuildCurve();
    emit pointsChanged(m_points.size());
}

// ---------------------------------------------------------------------------
// 顶点切线编辑
// ---------------------------------------------------------------------------

curve::VertexTangent::Mode PlotWidget::selectedVertexMode() const {
    if (m_selectedIndex >= 0 &&
        m_selectedIndex < static_cast<int>(m_tangentControls.size()))
        return m_tangentControls[static_cast<std::size_t>(m_selectedIndex)].mode;
    return m_pendingVertexMode;
}

void PlotWidget::setSelectedVertexMode(curve::VertexTangent::Mode mode) {
    m_pendingVertexMode = mode;  // 记住，供下次右键选新节点时套用
    if (m_selectedIndex >= 0 && m_selectedIndex < static_cast<int>(m_points.size()))
        applyVertexModeTo(m_selectedIndex, mode);
}

void PlotWidget::selectNode(int index) {
    m_selectedIndex = index;
    if (index >= 0 && index < static_cast<int>(m_points.size())) {
        auto& c = m_tangentControls[static_cast<std::size_t>(index)];
        // 新节点（未编辑）选中即套用当前"顶点模式"
        if (c.mode == curve::VertexTangent::Mode::Free &&
            m_pendingVertexMode != curve::VertexTangent::Mode::Free)
            applyVertexModeTo(index, m_pendingVertexMode);
    }
    emit nodeSelected(index);
    update();
}

void PlotWidget::applyVertexModeTo(int index, curve::VertexTangent::Mode mode) {
    auto& c = m_tangentControls[static_cast<std::size_t>(index)];
    if (mode == curve::VertexTangent::Mode::Free) {
        c = curve::VertexTangent{};  // 重置为未编辑
    } else if (mode == curve::VertexTangent::Mode::Smooth) {
        if (c.mode != curve::VertexTangent::Mode::Smooth) {
            if (c.mode == curve::VertexTangent::Mode::Straight ||
                c.mode == curve::VertexTangent::Mode::Corner)
                c.t = c.tRight;            // 从分裂切线进入：沿用右侧
            else
                c.t = estimateTangent(index);  // 从未编辑进入：用当前曲线切线估计
            c.mode = curve::VertexTangent::Mode::Smooth;
        }
    } else {  // Straight / Corner
        if (c.mode == curve::VertexTangent::Mode::Smooth) {
            c.tLeft = c.t;
            c.tRight = c.t;
        } else if (c.mode == curve::VertexTangent::Mode::Free) {
            const curve::Point est = estimateTangent(index);
            c.tLeft = est;
            c.tRight = est;
        }
        c.mode = mode;
        if (mode == curve::VertexTangent::Mode::Straight)
            straightenLeft(c);  // 进入直线模式时强制两侧共线（取右柄方向为直线方向）
    }
    rebuildCurve();
}

// 用弦长参数化的中央差分估计节点 k 的切线 dP/dt（与 parameterizeChordal 归一化一致）
curve::Point PlotWidget::estimateTangent(int index) const {
    const int n = static_cast<int>(m_points.size());
    if (n < 2) return {};
    std::vector<double> t(static_cast<std::size_t>(n), 0.0);
    double acc = 0.0;
    for (int i = 1; i < n; ++i) {
        const QPointF d = m_points[i] - m_points[i - 1];
        acc += std::hypot(d.x(), d.y());
        t[static_cast<std::size_t>(i)] = acc;
    }
    if (acc <= 0.0) return {};
    for (double& tv : t) tv /= acc;

    const auto slope = [&](int a, int b) -> QPointF {
        const double dt = t[static_cast<std::size_t>(b)] - t[static_cast<std::size_t>(a)];
        if (dt <= 0.0) return QPointF();
        const QPointF dp = m_points[b] - m_points[a];
        return QPointF(dp.x() / dt, dp.y() / dt);
    };
    QPointF m;
    if (index == 0)           m = slope(0, 1);
    else if (index == n - 1)  m = slope(n - 2, n - 1);
    else                      m = slope(index - 1, index + 1);
    return {m.x(), m.y()};
}

bool PlotWidget::handleEndpoints(QPointF& rightEnd, QPointF& leftEnd) const {
    if (m_selectedIndex < 0 || m_selectedIndex >= static_cast<int>(m_points.size()))
        return false;
    const curve::VertexTangent& c =
        m_tangentControls[static_cast<std::size_t>(m_selectedIndex)];
    if (c.mode == curve::VertexTangent::Mode::Free) return false;
    const QPointF pk = dataToWidget(m_points[m_selectedIndex]);
    if (c.mode == curve::VertexTangent::Mode::Smooth) {
        const QPointF tv = dataVecToWidget(
            QPointF(c.t.x * kHandleDisplayScale, c.t.y * kHandleDisplayScale));
        rightEnd = pk + tv;
        leftEnd = pk - tv;
    } else {
        rightEnd = pk + dataVecToWidget(
            QPointF(c.tRight.x * kHandleDisplayScale, c.tRight.y * kHandleDisplayScale));
        leftEnd = pk - dataVecToWidget(
            QPointF(c.tLeft.x * kHandleDisplayScale, c.tLeft.y * kHandleDisplayScale));
    }
    return true;
}

QPointF PlotWidget::dataVecToWidget(const QPointF& v) const {
    const QRectF r = dataRect();
    return QPointF(v.x() * r.width(), -v.y() * r.height());  // y 轴翻转
}

// 直线顶点：拖动右柄后，左柄切线保持与右柄同向共线（G¹ 要求进/出切线同向），长度不变
void PlotWidget::straightenLeft(curve::VertexTangent& c) {
    const double lenR = std::hypot(c.tRight.x, c.tRight.y);
    const double lenL = std::hypot(c.tLeft.x, c.tLeft.y);
    if (lenR < 1e-9) return;
    const double s = lenL / lenR;
    c.tLeft = {s * c.tRight.x, s * c.tRight.y};
}

// 直线顶点：拖动左柄后，右柄切线保持与左柄同向共线，长度不变
void PlotWidget::straightenRight(curve::VertexTangent& c) {
    const double lenL = std::hypot(c.tLeft.x, c.tLeft.y);
    const double lenR = std::hypot(c.tRight.x, c.tRight.y);
    if (lenL < 1e-9) return;
    const double s = lenR / lenL;
    c.tRight = {s * c.tLeft.x, s * c.tLeft.y};
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

    // 分段 Bezier 的中间控制点（仅当该层勾选时显示）：
    // 浅色虚线控制多边形 P_i — B1 — B2 — P_{i+1} + 靛蓝方块（可左键拖动）
    if (isLayerVisible(QStringLiteral("bezier_catmull_rom")) &&
        m_bezierSegments.size() + 1 == static_cast<std::size_t>(m_points.size())) {
        painter.setPen(QPen(QColor(0x4b, 0x00, 0x82, 110), 1, Qt::DashLine));
        for (std::size_t s = 0; s < m_bezierSegments.size(); ++s) {
            const QPointF p0 = dataToWidget(m_points[s]);
            const QPointF p1 = dataToWidget(m_points[s + 1]);
            const QPointF b1 = bezierControlWidgetPos(static_cast<int>(2 * s));
            const QPointF b2 = bezierControlWidgetPos(static_cast<int>(2 * s + 1));
            QPolygonF poly;
            poly << p0 << b1 << b2 << p1;
            painter.drawPolyline(poly);
        }
        painter.setPen(QPen(Qt::black, 1));
        painter.setBrush(QColor(0x4b, 0x00, 0x82));
        for (std::size_t s = 0; s < m_bezierSegments.size(); ++s) {
            for (int k = 0; k < 2; ++k) {
                const QPointF w = bezierControlWidgetPos(static_cast<int>(2 * s + k));
                painter.drawRect(QRectF(w.x() - 4, w.y() - 4, 8, 8));
            }
        }
        // 拖动中的控制点：外圈描边高亮
        if (m_dragBezier >= 0) {
            const QPointF w = bezierControlWidgetPos(m_dragBezier);
            if (!w.isNull()) {
                painter.setPen(QPen(Qt::black, 2));
                painter.setBrush(Qt::NoBrush);
                painter.drawRect(QRectF(w.x() - 6, w.y() - 6, 12, 12));
            }
        }
    }

    // 数据点（红）
    painter.setPen(QPen(Qt::black, 1));
    painter.setBrush(Qt::red);
    for (const QPointF& p : m_points) {
        painter.drawEllipse(dataToWidget(p), 4.0, 4.0);
    }

    // 拖动中的点：外圈描边高亮（选中标志）
    if (m_dragIndex >= 0 && m_dragIndex < m_points.size()) {
        painter.setPen(QPen(QColor(0x1f, 0x77, 0xb4), 2));  // 蓝描边
        painter.setBrush(Qt::NoBrush);
        painter.drawEllipse(dataToWidget(m_points[m_dragIndex]), 7.0, 7.0);
    }

    // 选中的节点（右键）：绿色外圈 + 切线柄
    if (m_selectedIndex >= 0 && m_selectedIndex < static_cast<int>(m_points.size())) {
        const QPointF pk = dataToWidget(m_points[m_selectedIndex]);
        painter.setPen(QPen(QColor(0x2e, 0x8b, 0x57), 2));  // 绿描边
        painter.setBrush(Qt::NoBrush);
        painter.drawEllipse(pk, 7.0, 7.0);

        QPointF rightEnd, leftEnd;
        if (handleEndpoints(rightEnd, leftEnd)) {
            painter.setPen(QPen(QColor(0x1f, 0x77, 0xb4), 1.5));  // 右柄（蓝）
            painter.drawLine(pk, rightEnd);
            painter.setPen(QPen(QColor(0xff, 0x7f, 0x0e), 1.5));  // 左柄（橙）
            painter.drawLine(pk, leftEnd);
            painter.setPen(Qt::black);
            painter.setBrush(QColor(0x1f, 0x77, 0xb4));
            painter.drawEllipse(rightEnd, 3.0, 3.0);
            painter.setBrush(QColor(0xff, 0x7f, 0x0e));
            painter.drawEllipse(leftEnd, 3.0, 3.0);
        }
    }

    // 顶部提示 + 图例
    QFont f = font();
    f.setPointSize(10);
    painter.setFont(f);
    const QRectF hint(dr.left(), 6, dr.width(), 20);
    painter.setPen(QColor(0x88, 0x88, 0x88));
    painter.drawText(hint, Qt::AlignLeft | Qt::AlignVCenter,
                     QStringLiteral("左键加点 · 左键长按拖点 · 右键选点出切线柄 · "
                                    "左键拖靛蓝方块=Bezier控制点 · 当前 %1 个点")
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
