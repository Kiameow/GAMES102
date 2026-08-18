#include <QApplication>
#include <QString>

#include "MainWindow.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    app.setApplicationName(QStringLiteral("GAMES102 Curve Lab"));
    app.setOrganizationName(QStringLiteral("GAMES102"));

    MainWindow window;
    window.show();
    return app.exec();
}
