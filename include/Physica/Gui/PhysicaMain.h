/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_PHYSICAMAIN_H
#define PHYSICA_PHYSICAMAIN_H

#include <QtWidgets/QMainWindow>
#include "Physica/Gui/Editor/EditorMain.h"
#include "Physica/Gui/Menu/Calculator.h"

namespace Physica::Gui {
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
}


#endif