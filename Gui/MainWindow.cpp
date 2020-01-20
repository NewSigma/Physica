#include "Header/MainWindow.h"
#include "Header/CubePlayer.h"
#include <QGuiApplication>
#include <QScreen>

MainWindow::MainWindow() {
    QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
    resize(primaryScreenRec.width(), primaryScreenRec.height());
    setWindowTitle("Physica");
    icon = new QIcon("../Resources/icon.png");
    setWindowIcon(*icon);

    format = new QSurfaceFormat();
    format->setSamples(4);
    QSurfaceFormat::setDefaultFormat(*format);

    plotter = new CubePlayer(this);
    plotter->reloadTimer(20);
    setCentralWidget(plotter);
}

MainWindow::~MainWindow() {
    delete icon;
    delete format;
    destroy();
}