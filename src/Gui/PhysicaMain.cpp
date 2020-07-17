/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QApplication>
#include "Physica/Gui/PhysicaMain.h"
#include "QMenuBar"
#include "QScreen"
#include "QKeyEvent"
#include "Physica/Gui/Settings.h"
#include "Physica/Gui/LongtonAnt.h"
#include "Physica/Gui/GlobalStyle.h"

namespace Physica::Gui {
    PhysicaMain::PhysicaMain() : lastCalcTime(0), calculator(nullptr) {
        /* Resize */ {
            QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
            resize(primaryScreenRec.width(), primaryScreenRec.height());
        }
        //Basic settings
        setWindowTitle(QApplication::applicationName() + ' ' + QApplication::applicationVersion() + "[*]");
        QIcon icon = QIcon(":/icon.png");
        setAttribute(Qt::WA_DeleteOnClose);
        setWindowIcon(icon);
        /* Menu */ {
            menuBar = new QMenuBar(this);
            setMenuBar(menuBar);
            file = menuBar->addMenu("File");
            file_new = file->addAction("New");
            file_open = file->addAction("Open");
            file->addSeparator();
            file_settings = file->addAction("Settings");
            connect(file_settings, &QAction::triggered, [&]() { new Settings(this); });

            insert = menuBar->addMenu("Insert");
            insert_cell = insert->addAction("Insert Cell");
            connect(insert_cell, &QAction::triggered, [&]() { textEdit->textCursor().insertFrame(QTextFrameFormat()); });

            view = menuBar->addMenu("View");
            view_fullscreen = view->addAction("Full Screen");
            connect(view_fullscreen, SIGNAL(triggered()), SLOT(on_click_fullscreen()));

            evaluate = menuBar->addMenu("Evaluate");
            evaluate_calculator = evaluate->addAction("Calculator");
            connect(evaluate_calculator, SIGNAL(triggered()), SLOT(on_click_calculator()));
            //InT stands by Interesting Things.
            InT = menuBar->addMenu("InT");
            lontonAnt = InT->addAction("Lonton Ant");
            connect(lontonAnt, &QAction::triggered, []() { new LongtonAnt(); });
        }
        //Central
        textEdit = new EditorMain(this);
        textEdit->installEventFilter(this);
        connect(textEdit->document(), SIGNAL(contentsChanged()), SLOT(on_modified()));
        setCentralWidget(textEdit);
        //StatusBar
        statusBar()->showMessage("Calculate Time: 0 ms");

        show();
    }
////////////////////////////////////////////Events////////////////////////////////////////////
    bool PhysicaMain::eventFilter(QObject* obj, QEvent* event) {
        if (obj == textEdit) {
            if (event->type() == QEvent::KeyPress) {
                auto keyEvent = (QKeyEvent*)event;
                if(keyEvent->modifiers() == Qt::ShiftModifier && keyEvent->key() == Qt::Key_Return) {
                    lastCalcTime = clock();
                    //TODO
                    clock_t end = clock();
                    statusBar()->showMessage(QString::fromStdString("Calculate Time: " + std::to_string(end - lastCalcTime) + " ms"));
                    return true;
                }
            }
            return false;
        }
        else
            return QObject::eventFilter(obj, event);
    }

    void PhysicaMain::keyPressEvent(QKeyEvent* event) {
        if(event->modifiers() == Qt::NoModifier) {
            switch(event->key()) {
                case Qt::Key_F11:
                    on_click_fullscreen();
                default:;
            }
        }
        else {
            switch(event->modifiers()) {
                default:;
            }
        }
    }
////////////////////////////////////////////SLOTS////////////////////////////////////////////
    void PhysicaMain::on_modified() {
        setWindowModified(true);
        disconnect(textEdit->document(), &QTextDocument::contentsChanged, this, &PhysicaMain::on_modified);
    }

    void PhysicaMain::on_click_fullscreen() {
        if(windowState() == Qt::WindowFullScreen)
            setWindowState(Qt::WindowMaximized);
        else
            setWindowState(Qt::WindowFullScreen);
    }

    void PhysicaMain::on_click_calculator() {
        if(calculator == nullptr) {
            calculator = new Calculator();
            connect(calculator, &Calculator::destroyed, [&]() { calculator = nullptr; });
        }
    }
}