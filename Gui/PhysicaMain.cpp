#include <QApplication>
#include "Header/PhysicaMain.h"
#include "QMenuBar"
#include "QScreen"
#include "QKeyEvent"
#include "Header/Settings.h"
#include "Header/LongtonAnt.h"
#include "Header/GlobalStyle.h"

/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
PhysicaMain::PhysicaMain() : lastCalcTime(0), calculator(nullptr) {
    /* Resize */ {
        QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
        resize(primaryScreenRec.width(), primaryScreenRec.height());
    }
    //Basic settings
    setWindowTitle(QApplication::applicationName() + "[*]");
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
        connect(file_settings, SIGNAL(triggered()), SLOT(on_click_settings()));

        view = menuBar->addMenu("View");
        view_fullscreen = view->addAction("Full Screen");
        connect(view_fullscreen, SIGNAL(triggered()), SLOT(on_click_fullscreen()));

        evaluate = menuBar->addMenu("Evaluate");
        evaluate_calculator = evaluate->addAction("Calculator");
        connect(evaluate_calculator, SIGNAL(triggered()), SLOT(on_click_calculator()));
        //InT stands by Interesting Things.
        InT = menuBar->addMenu("InT");
        lontonAnt = InT->addAction("Lonton Ant");
        connect(lontonAnt, SIGNAL(triggered()), SLOT(on_click_lontonAnt()));
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
                statusBar()->showMessage(QString::fromStdString("Calculate Time: " + std::to_string(clock() - lastCalcTime) + " ms"));
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
    disconnect(textEdit->document(), SIGNAL(contentsChanged()), this, SLOT(on_modified()));
}

void PhysicaMain::on_click_settings() {
    new Settings(this);
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
        connect(calculator, SIGNAL(destroyed()), SLOT(on_calculator_closed()));
    }
}

void PhysicaMain::on_calculator_closed() {
    calculator = nullptr;
}
////////////////////////////////////////////InT////////////////////////////////////////////
void PhysicaMain::on_click_lontonAnt() {
    new LongtonAnt();
}