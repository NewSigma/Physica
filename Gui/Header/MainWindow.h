#ifndef PHYSICA_C_MAINWINDOW_H
#define PHYSICA_C_MAINWINDOW_H

#include <QMainWindow>
#include <QtGui/QSurfaceFormat>
#include <QtWidgets/QLayout>
#include "Plotter.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow();
    ~MainWindow();
private:
    QIcon* icon;
    QSurfaceFormat* format;
    Plotter* plotter;
};

#endif