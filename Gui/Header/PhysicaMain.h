#ifndef PHYSICA_PHYSICAMAIN_H
#define PHYSICA_PHYSICAMAIN_H

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QTextEdit>
#include "Calculator.h"

class PhysicaMain : public QMainWindow {
    Q_OBJECT
public:
    PhysicaMain();
    ~PhysicaMain() override;
private:
    QIcon* icon;
    QMenuBar* menuBar;
        QMenu* files;
            QAction* files_new;
            QAction* files_open;
            QAction* files_settings;
        QMenu* tools;
            QAction* tools_calculator;
    QTextEdit* textEdit;
    Calculator* calculator;
private slots:
    void on_clicked_calculator();
    void on_calculator_closed();
};

#endif