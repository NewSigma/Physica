#include <QtGui/QGuiApplication>
#include "Header/PhysicaMain.h"
#include "QMenuBar"
#include "QScreen"

PhysicaMain::PhysicaMain() : calculator(nullptr) {
    /* Resize */ {
        QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
        resize(primaryScreenRec.width() * 3 / 4, primaryScreenRec.height() * 3 / 4);
    }
    //Basic settings
    setWindowTitle("Physica");
    icon = new QIcon(":/icon.png");
    setWindowIcon(*icon);
    /* Menu */ {
        menuBar = new QMenuBar();
        setMenuBar(menuBar);
        files = menuBar->addMenu(tr("Files"));
        files_new = files->addAction("New");
        files_open = files->addAction("Open");
        files->addSeparator();
        files_settings = files->addAction(tr("Settings"));
        tools = menuBar->addMenu(tr("Tools"));
        tools_calculator = tools->addAction(tr("Calculator"));
        connect(tools_calculator, SIGNAL(triggered()), SLOT(on_clicked_calculator()));
    }
    textEdit = new QTextEdit();
    setCentralWidget(textEdit);

    show();
}

PhysicaMain::~PhysicaMain() {
    delete icon;
    delete menuBar;
    delete textEdit;
}

void PhysicaMain::on_clicked_calculator() {
    if(calculator == nullptr) {
        calculator = new Calculator();
        connect(calculator, SIGNAL(destroyed()), SLOT(on_calculator_closed()));
    }
}

void PhysicaMain::on_calculator_closed() {
    calculator = nullptr;
}