#ifndef PHYSICA_MAINWINDOW_H
#define PHYSICA_MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include "Calculator.h"

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();
    ~MainWindow() override;
private:
    QIcon* icon;
    QPushButton* button;
    Calculator* calculator;
private slots:
    void on_clicked();
};

#endif