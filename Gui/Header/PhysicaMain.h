#ifndef PHYSICA_PHYSICAMAIN_H
#define PHYSICA_PHYSICAMAIN_H

#include <QtWidgets/QMainWindow>
#include "EditorMain.h"
#include "Calculator.h"

class PhysicaMain : public QMainWindow {
    Q_OBJECT
private:
    long lastCalcTime;

    QMenuBar* menuBar;
        QMenu* file;
            QAction* file_new;
            QAction* file_open;
            QAction* file_settings;
        QMenu* insert;
            QAction* insert_cell;
        QMenu* view;
            QAction* view_fullscreen;
        QMenu* evaluate;
            QAction* evaluate_calculator;
        QMenu* InT;
            QAction* lontonAnt;
    QTextEdit* textEdit;
    Calculator* calculator;
public:
    PhysicaMain();
private:
    ////////////////////////////////////////////Events////////////////////////////////////////////
    bool eventFilter(QObject* obj, QEvent* event) override;
    void keyPressEvent(QKeyEvent* event) override;
private slots:
    ////////////////////////////////////////////SLOTS////////////////////////////////////////////
    void on_modified();
    void on_click_fullscreen();
    void on_click_calculator();
};

#endif