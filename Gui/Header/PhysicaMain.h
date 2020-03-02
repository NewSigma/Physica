#ifndef PHYSICA_PHYSICAMAIN_H
#define PHYSICA_PHYSICAMAIN_H

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QTextEdit>
#include "Calculator.h"

class PhysicaMain : public QMainWindow {
    Q_OBJECT
private:
    QIcon* icon;
    QMenuBar* menuBar;
        QMenu* file;
            QAction* file_new;
            QAction* file_open;
            QAction* file_settings;
        QMenu* view;
            QAction* view_fullscreen;
        QMenu* evaluate;
            QAction* evaluate_calculator;
    QTextEdit* textEdit;
    Calculator* calculator;
public:
    PhysicaMain();
    ~PhysicaMain() override;
private:
    bool eventFilter(QObject* obj, QEvent* event);
    void keyPressEvent(QKeyEvent* event) override;
private slots:
    void on_click_fullscreen();
    void on_click_calculator();
    void on_calculator_closed();
};

#endif