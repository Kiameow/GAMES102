#include "MainWindow.h"

#include "math/CurveMath.h"

#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QStatusBar>
#include <QVBoxLayout>

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle(QStringLiteral("GAMES102 — 曲线插值与逼近"));
    resize(1360, 860);

    m_plot = new PlotWidget(this);

    // ---- 控制面板 ----
    m_spinDegree = new QSpinBox(this);
    m_spinDegree->setRange(1, 20);
    m_spinDegree->setValue(3);
    m_spinDegree->setToolTip(QStringLiteral("最小二乘逼近多项式的次数"));

    m_spinSigma = new QDoubleSpinBox(this);
    m_spinSigma->setRange(0.05, 2.0);
    m_spinSigma->setDecimals(2);
    m_spinSigma->setSingleStep(0.05);
    m_spinSigma->setValue(0.3);
    m_spinSigma->setToolTip(QStringLiteral("高斯基函数的宽度 σ（越小曲线越局部）"));

    m_spinLambda = new QDoubleSpinBox(this);
    m_spinLambda->setRange(0.0, 5.0);
    m_spinLambda->setDecimals(2);
    m_spinLambda->setSingleStep(0.05);
    m_spinLambda->setValue(0.1);
    m_spinLambda->setToolTip(QStringLiteral("岭回归正则系数 λ（越大系数越被压向 0、曲线越平滑；0 即普通最小二乘）"));

    auto* btnShowAll = new QPushButton(QStringLiteral("显示所有"), this);
    auto* btnHideAll = new QPushButton(QStringLiteral("隐藏全部"), this);

    auto* btnClear = new QPushButton(QStringLiteral("清除所有点"), this);
    btnClear->setStyleSheet(QStringLiteral("color:#c0392b;"));

    auto* degreeRow = new QHBoxLayout;
    degreeRow->addWidget(new QLabel(QStringLiteral("逼近次数:"), this));
    degreeRow->addWidget(m_spinDegree, 1);

    auto* sigmaRow = new QHBoxLayout;
    sigmaRow->addWidget(new QLabel(QStringLiteral("高斯σ:"), this));
    sigmaRow->addWidget(m_spinSigma, 1);

    auto* lambdaRow = new QHBoxLayout;
    lambdaRow->addWidget(new QLabel(QStringLiteral("岭回归λ:"), this));
    lambdaRow->addWidget(m_spinLambda, 1);

    auto* info = new QLabel(
        QStringLiteral("左键：添加红点\n右键：删除最近的红点\n\n"
                       "勾选按钮可叠加显示多条曲线，\n"
                       "「显示所有」一次全开，不同颜色\n区分不同算法。"),
        this);
    info->setWordWrap(true);
    info->setStyleSheet(QStringLiteral("color:#555;"));

    auto* controls = new QVBoxLayout;
    controls->addWidget(new QLabel(QStringLiteral("<b>GAMES102 作业工具</b>"), this));
    controls->addSpacing(8);
    controls->addLayout(degreeRow);
    controls->addLayout(sigmaRow);
    controls->addLayout(lambdaRow);

    // ---- 注册算法（以后新增算法只需在这里加一行 addAlgorithm）----
    // 曲线颜色参考（高区分度备选色板，避免相邻算法撞色）：
    //   蓝 #1f77b4   橙 #ff7f0e   绿 #2ca02c   紫 #9467bd
    //   青 #17becf   棕 #8c564b   粉 #e377c2   灰 #7f7f7f
    // 注意：红 #d62728 与数据点(红色)冲突，曲线颜色慎用红色系。
    addAlgorithm(QStringLiteral("lagrange"), QStringLiteral("插值-Lagrange"),
                 QColor(0x1f, 0x77, 0xb4),
                 [](const std::vector<curve::Point>& pts, int samples) {
                     return curve::lagrangeInterpolate(pts, samples);
                 },
                 controls);

    addAlgorithm(QStringLiteral("power"), QStringLiteral("插值-幂基"),
                 QColor(0xff, 0x7f, 0x0e),
                 [](const std::vector<curve::Point>& pts, int samples) {
                     return curve::powerBasePolynomialInterpolate(pts, samples);
                 },
                 controls);

    addAlgorithm(QStringLiteral("gauss"), QStringLiteral("插值-Gauss基"),
                 QColor(0x94, 0x67, 0xbd),
                 [this](const std::vector<curve::Point>& pts, int samples) {
                     return curve::gaussBasePolynomialInterpolate(pts, m_plot->gaussianSigma(),
                                                                  samples);
                 },
                 controls);

    addAlgorithm(QStringLiteral("least_squares"), QStringLiteral("逼近-幂函数最小二乘"),
                 QColor(0x2c, 0xa0, 0x2c),
                 [this](const std::vector<curve::Point>& pts, int samples) {
                     return curve::leastSquaresPolynomial(pts, m_plot->approximationDegree(),
                                                          samples);
                 },
                 controls);

    addAlgorithm(QStringLiteral("ridge"), QStringLiteral("逼近-岭回归"),
                 QColor(0x17, 0xbe, 0xcf),
                 [this](const std::vector<curve::Point>& pts, int samples) {
                     return curve::ridgeRegression(pts, m_plot->approximationDegree(),
                                                   m_plot->lambda(), samples);
                 },
                 controls);

    auto* allRow = new QHBoxLayout;
    allRow->addWidget(btnShowAll);
    allRow->addWidget(btnHideAll);
    controls->addLayout(allRow);
    controls->addWidget(btnClear);
    controls->addSpacing(12);
    controls->addWidget(info);
    controls->addStretch(1);

    auto* panel = new QWidget(this);
    panel->setLayout(controls);
    panel->setFixedWidth(280);

    auto* layout = new QHBoxLayout;
    layout->addWidget(m_plot, 1);
    layout->addWidget(panel);
    auto* central = new QWidget(this);
    central->setLayout(layout);
    setCentralWidget(central);

    // ---- 信号连接 ----
    connect(btnShowAll, &QPushButton::clicked, this, &MainWindow::onShowAllClicked);
    connect(btnHideAll, &QPushButton::clicked, this, &MainWindow::onHideAllClicked);
    connect(btnClear, &QPushButton::clicked, this, &MainWindow::onClearClicked);
    connect(m_plot, &PlotWidget::pointsChanged, this, &MainWindow::onPointsChanged);
    connect(m_spinDegree, QOverload<int>::of(&QSpinBox::valueChanged), this,
            &MainWindow::onDegreeChanged);
    connect(m_spinSigma, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &MainWindow::onSigmaChanged);
    connect(m_spinLambda, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &MainWindow::onLambdaChanged);

    onPointsChanged(0);
}

// 通用注册入口：注册算法层 + 生成按钮，所有冗余逻辑都在这里，加新算法不用再碰
void MainWindow::addAlgorithm(const QString& key, const QString& name, const QColor& color,
                              CurveCompute compute, QVBoxLayout* layout) {
    // 1. 注册到画布（名字 / 颜色 / 计算函数）
    CurveLayer layer;
    layer.key = key;
    layer.name = name;
    layer.color = color;
    layer.compute = std::move(compute);
    m_plot->registerLayer(layer);

    // 2. 生成一个可勾选按钮，勾选 = 显示该算法曲线
    auto* btn = new QPushButton(name, this);
    btn->setCheckable(true);
    layout->addWidget(btn);
    connect(btn, &QPushButton::toggled, this, [this, key, name](bool on) {
        m_plot->setLayerVisible(key, on);
        statusBar()->showMessage(on ? QStringLiteral("已显示：%1").arg(name)
                                    : QStringLiteral("已隐藏：%1").arg(name),
                                 3000);
    });
    m_algorithmButtons.insert(key, btn);
}

void MainWindow::onShowAllClicked() {
    for (QPushButton* btn : m_algorithmButtons) btn->setChecked(true);
    statusBar()->showMessage(QStringLiteral("已显示所有算法结果"), 2000);
}

void MainWindow::onHideAllClicked() {
    for (QPushButton* btn : m_algorithmButtons) btn->setChecked(false);
    statusBar()->showMessage(QStringLiteral("已隐藏所有算法结果"), 2000);
}

void MainWindow::onClearClicked() {
    m_plot->clearPoints();
    for (QPushButton* btn : m_algorithmButtons) btn->setChecked(false);
    statusBar()->showMessage(QStringLiteral("已清除所有点"), 2000);
}

void MainWindow::onPointsChanged(int count) {
    statusBar()->showMessage(
        QStringLiteral("当前 %1 个数据点（至少 2 个点才能显示曲线）").arg(count), 3000);
}

void MainWindow::onDegreeChanged(int degree) {
    m_plot->setApproximationDegree(degree);
    if (m_plot->isLayerVisible(QStringLiteral("least_squares"))) {
        statusBar()->showMessage(QStringLiteral("逼近次数已调整为 %1").arg(degree), 2000);
    }
}

void MainWindow::onSigmaChanged(double sigma) {
    m_plot->setGaussianSigma(sigma);
    if (m_plot->isLayerVisible(QStringLiteral("gauss"))) {
        statusBar()->showMessage(QStringLiteral("高斯σ已调整为 %1").arg(sigma), 2000);
    }
}

void MainWindow::onLambdaChanged(double lambda) {
    m_plot->setLambda(lambda);
    if (m_plot->isLayerVisible(QStringLiteral("ridge"))) {
        statusBar()->showMessage(QStringLiteral("岭回归λ已调整为 %1").arg(lambda), 2000);
    }
}
