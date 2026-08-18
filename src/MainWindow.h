#pragma once

#include <QMainWindow>
#include <QMap>
#include <QString>

#include "PlotWidget.h"

class QPushButton;
class QDoubleSpinBox;
class QSpinBox;
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
    void onLambdaChanged(double lambda);
    void onShowAllClicked();
    void onHideAllClicked();

private:
    // 通用注册入口：往画布注册算法层 + 生成一个可勾选按钮。
    // 以后新增算法只需在构造函数里调用一次 addAlgorithm(...)。
    void addAlgorithm(const QString& key, const QString& name, const QColor& color,
                      CurveCompute compute, QVBoxLayout* layout);

    PlotWidget* m_plot = nullptr;
    QMap<QString, QPushButton*> m_algorithmButtons;  // 算法 key -> 按钮
    QSpinBox* m_spinDegree = nullptr;
    QDoubleSpinBox* m_spinSigma = nullptr;
    QDoubleSpinBox* m_spinLambda = nullptr;
};
