#include <QApplication>
#include <iostream>
#include "Header/PhysicaMain.h"
#include "QMenuBar"
#include "QScreen"
#include "QKeyEvent"

PhysicaMain::PhysicaMain() : calculator(nullptr) {
    /* Resize */ {
        QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
        resize(primaryScreenRec.width(), primaryScreenRec.height());
    }
    //Basic settings
    setWindowTitle(QApplication::applicationName() + " " + QApplication::applicationVersion());
    icon = new QIcon(":/icon.png");
    setWindowIcon(*icon);
    /* Menu */ {
        menuBar = new QMenuBar();
        setMenuBar(menuBar);
        file = menuBar->addMenu("File");
        file_new = file->addAction("New");
        file_open = file->addAction("Open");
        file->addSeparator();
        file_settings = file->addAction("Settings");

        view = menuBar->addMenu("View");
        view_fullscreen = view->addAction("Full Screen");
        connect(view_fullscreen, SIGNAL(triggered()), SLOT(on_click_fullscreen()));

        evaluate = menuBar->addMenu("Evaluate");
        evaluate_calculator = evaluate->addAction("Calculator");
        connect(evaluate_calculator, SIGNAL(triggered()), SLOT(on_click_calculator()));
    }
    textEdit = new QTextEdit();
    textEdit->installEventFilter(this);
    setCentralWidget(textEdit);

    show();
}

PhysicaMain::~PhysicaMain() {
    delete icon;
    delete menuBar;
    delete textEdit;
}

bool PhysicaMain::eventFilter(QObject* obj, QEvent* event) {
    if (obj == textEdit) {
        if (event->type() == QEvent::KeyPress) {
            auto keyEvent = (QKeyEvent*)event;
            if(keyEvent->modifiers() == Qt::ShiftModifier && keyEvent->key() == Qt::Key_Return) {
                textEdit->insertPlainText(tr("Done"));
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