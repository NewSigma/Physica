#include "Header/MainWindow.h"

MainWindow::MainWindow() {
    //Basic settings
    setWindowTitle("Physica");
    resize(1440, 810);
    icon = new QIcon(":/icon.png");
    setWindowIcon(*icon);

    button = new QPushButton("Calculator");
    setCentralWidget(button);
    calculator = nullptr;
    connect(button, SIGNAL(clicked()), SLOT(on_clicked()));

    show();
}

MainWindow::~MainWindow() {
    delete icon;
    delete button;
    delete calculator;
}

void MainWindow::on_clicked() {
    if(!calculator)
        calculator = new Calculator();
}