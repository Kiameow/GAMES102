#pragma once

#include <QJsonObject>
#include <QMainWindow>
#include <QMap>
#include <QProcess>
#include <QString>

#include "PlotWidget.h"

class QComboBox;
class QDoubleSpinBox;
class QPushButton;
class QSpinBox;
class QTimer;
class QVBoxLayout;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);

private slots:
    void onClearClicked();
    void onPointsChanged(int count);
    void onDegreeChanged(int degree);
    void onSigmaChanged(double sigma);
    void onGaussCentersChanged(int centers);
    void onParamBasisChanged(int index);
    void onVertexModeChanged(int index);
    void onNodeSelected(int index);
    void onLambdaChanged(double lambda);
    void onShowAllClicked();
    void onHideAllClicked();
    void onRbfToggled(bool on);
    void onRbfFinished(int exitCode, QProcess::ExitStatus status);
    void onCaptureScreenshot();

private:
    // 通用注册入口：往画布注册算法层 + 生成一个可勾选按钮。
    // 以后新增算法只需在构造函数里调用一次 addAlgorithm(...)。
    void addAlgorithm(const QString& key, const QString& name, const QColor& color,
                      CurveCompute compute, QVBoxLayout* layout);
    // 按当前"参数拟合基函数"选择调用对应拟合（幂基最小二乘 / Gauss 基最小二乘）
    std::vector<curve::Point> fitParametricSelected(const std::vector<curve::Point>& pts,
                                                    const std::vector<double>& t, int samples);
    // RBF：找项目根目录（向上找 CMakeLists.txt）
    QString findProjectRoot() const;
    void startRbfTraining();
    bool readRbfOutput(QJsonObject& obj) const;  // 读 build/io/output.json，失败返回 false
    void onRbfPoll();                            // 训练中定时轮询，实时刷新曲线

    PlotWidget* m_plot = nullptr;
    QMap<QString, QPushButton*> m_algorithmButtons;  // 算法 key -> 按钮
    QSpinBox* m_spinDegree = nullptr;
    QDoubleSpinBox* m_spinSigma = nullptr;
    QDoubleSpinBox* m_spinLambda = nullptr;
    QPushButton* m_btnRbf = nullptr;
    QPushButton* m_btnScreenshot = nullptr;
    QSpinBox* m_spinCenters = nullptr;
    QSpinBox* m_spinGaussCenters = nullptr;
    QSpinBox* m_spinEpochs = nullptr;
    QComboBox* m_comboParamBasis = nullptr;  // 参数拟合基函数：幂基 / Gauss基
    QComboBox* m_comboVertexMode = nullptr;  // 顶点模式：平滑 / 直线 / 角部 / 自由
    QProcess* m_rbfProcess = nullptr;
    QTimer* m_rbfTimer = nullptr;  // 训练中轮询 output.json 实现实时预览
};
