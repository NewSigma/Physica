#include "Header/Calculator.h"
#include "../Core/Header/ExprReader.h"
#include <QKeyEvent>
#include "../Core/Header/ElementaryFunction.h"

Calculator::Calculator() {
    //Const
    QFont font;
    font.setPointSize(20);
    //Main
    setAttribute(Qt::WA_DeleteOnClose);
    setWindowTitle("Calculator");
    setFixedSize(800,675);
    setWindowOpacity(0.9);
    setStyleSheet("background-color:rgb(50, 50, 50);");
    //Layout
    default_layout = new QVBoxLayout(this);
    //Editor
    editor_top = new QLineEdit();
    editor_top->setFixedSize(760, 46);
    editor_top->setStyleSheet("border:0px; color:rgb(255, 255, 255); background-color:rgb(50, 50, 50);");
    editor_top->setAlignment(Qt::AlignRight);
    editor_top->setReadOnly(true);
    editor_top->setFocusPolicy(Qt::NoFocus);
    default_layout->addWidget(editor_top);

    editor_bottom = new QLineEdit();
    editor_bottom->setFixedSize(760, 93);
    editor_bottom->setFont(font);
    editor_bottom->setStyleSheet("border:0px; color:rgb(255, 255, 255); background-color:rgb(50, 50, 50);");
    editor_bottom->setReadOnly(true);
    editor_bottom->setFocusPolicy(Qt::NoFocus);
    editor_bottom->setAlignment(Qt::AlignRight);
    default_layout->addWidget(editor_bottom);
    /* Buttons */ {
        layout = new QGridLayout();
        layout->setContentsMargins(2,3,2,3);
        default_layout->addLayout(layout);

        num0 = new QPushButton("0");
        num1 = new QPushButton("1");
        num2 = new QPushButton("2");
        num3 = new QPushButton("3");
        num4 = new QPushButton("4");
        num5 = new QPushButton("5");
        num6 = new QPushButton("6");
        num7 = new QPushButton("7");
        num8 = new QPushButton("8");
        num9 = new QPushButton("9");
        dot = new QPushButton(".");
        imagine = new QPushButton("i");
        add = new QPushButton("+");
        subtract = new QPushButton("-");
        multiply = new QPushButton("×");
        divide = new QPushButton("/");
        equal = new QPushButton("=");
        del = new QPushButton("Del");
        clear = new QPushButton("C");
        clear_entry = new QPushButton("CE");
        left_bracket = new QPushButton("(");
        right_bracket = new QPushButton(")");
        sqroot = new QPushButton("x²");
        square = new QPushButton("√x");
        percent = new QPushButton("%");

        num0->setFixedSize(150, 93);
        num1->setFixedSize(150, 93);
        num2->setFixedSize(150, 93);
        num3->setFixedSize(150, 93);
        num4->setFixedSize(150, 93);
        num5->setFixedSize(150, 93);
        num6->setFixedSize(150, 93);
        num7->setFixedSize(150, 93);
        num8->setFixedSize(150, 93);
        num9->setFixedSize(150, 93);
        dot->setFixedSize(150, 93);
        imagine->setFixedSize(150, 93);
        add->setFixedSize(150, 93);
        subtract->setFixedSize(150, 93);
        multiply->setFixedSize(150, 93);
        divide->setFixedSize(150, 93);
        equal->setFixedSize(150, 93);
        del->setFixedSize(150, 93);
        clear->setFixedSize(150, 93);
        clear_entry->setFixedSize(150, 93);
        left_bracket->setFixedSize(150, 93);
        right_bracket->setFixedSize(150, 93);
        sqroot->setFixedSize(150, 93);
        square->setFixedSize(150, 93);
        percent->setFixedSize(150, 93);

        auto sheet = "color:rgb(255, 255, 255); background-color:rgb(10, 10, 10);";
        auto sheet1 = "color:rgb(255, 255, 255); background-color:rgb(25, 25, 25);";
        num0->setStyleSheet(sheet);
        num1->setStyleSheet(sheet);
        num2->setStyleSheet(sheet);
        num3->setStyleSheet(sheet);
        num4->setStyleSheet(sheet);
        num5->setStyleSheet(sheet);
        num6->setStyleSheet(sheet);
        num7->setStyleSheet(sheet);
        num8->setStyleSheet(sheet);
        num9->setStyleSheet(sheet);
        dot->setStyleSheet(sheet1);
        imagine->setStyleSheet(sheet1);
        add->setStyleSheet(sheet1);
        subtract->setStyleSheet(sheet1);
        multiply->setStyleSheet(sheet1);
        divide->setStyleSheet(sheet1);
        equal->setStyleSheet(sheet1);
        del->setStyleSheet(sheet1);
        clear->setStyleSheet(sheet1);
        clear_entry->setStyleSheet(sheet1);
        left_bracket->setStyleSheet(sheet1);
        right_bracket->setStyleSheet(sheet1);
        sqroot->setStyleSheet(sheet1);
        square->setStyleSheet(sheet1);
        percent->setStyleSheet(sheet1);

        num0->setFont(font);
        num1->setFont(font);
        num2->setFont(font);
        num3->setFont(font);
        num4->setFont(font);
        num5->setFont(font);
        num6->setFont(font);
        num7->setFont(font);
        num8->setFont(font);
        num9->setFont(font);
        dot->setFont(font);
        imagine->setFont(font);
        add->setFont(font);
        subtract->setFont(font);
        multiply->setFont(font);
        divide->setFont(font);
        equal->setFont(font);
        del->setFont(font);
        clear->setFont(font);
        clear_entry->setFont(font);
        left_bracket->setFont(font);
        right_bracket->setFont(font);
        sqroot->setFont(font);
        square->setFont(font);
        percent->setFont(font);

        layout->addWidget(percent, 0, 0);
        layout->addWidget(clear_entry, 0, 1);
        layout->addWidget(clear, 0, 2);
        layout->addWidget(del, 0, 3);
        layout->addWidget(divide, 0, 4);
        layout->addWidget(square, 1, 0);
        layout->addWidget(num7, 1, 1);
        layout->addWidget(num8, 1, 2);
        layout->addWidget(num9, 1, 3);
        layout->addWidget(multiply, 1, 4);
        layout->addWidget(sqroot, 2, 0);
        layout->addWidget(num4, 2, 1);
        layout->addWidget(num5, 2, 2);
        layout->addWidget(num6, 2, 3);
        layout->addWidget(subtract, 2, 4);
        layout->addWidget(left_bracket, 3, 0);
        layout->addWidget(num1, 3, 1);
        layout->addWidget(num2, 3, 2);
        layout->addWidget(num3, 3, 3);
        layout->addWidget(add, 3, 4);
        layout->addWidget(right_bracket, 4, 0);
        layout->addWidget(imagine, 4, 1);
        layout->addWidget(num0, 4, 2);
        layout->addWidget(dot, 4, 3);
        layout->addWidget(equal, 4, 4);

        connect(num0, SIGNAL(clicked()), SLOT(on_click_num0()));
        connect(num1, SIGNAL(clicked()), SLOT(on_click_num1()));
        connect(num2, SIGNAL(clicked()), SLOT(on_click_num2()));
        connect(num3, SIGNAL(clicked()), SLOT(on_click_num3()));
        connect(num4, SIGNAL(clicked()), SLOT(on_click_num4()));
        connect(num5, SIGNAL(clicked()), SLOT(on_click_num5()));
        connect(num6, SIGNAL(clicked()), SLOT(on_click_num6()));
        connect(num7, SIGNAL(clicked()), SLOT(on_click_num7()));
        connect(num8, SIGNAL(clicked()), SLOT(on_click_num8()));
        connect(num9, SIGNAL(clicked()), SLOT(on_click_num9()));
        connect(dot, SIGNAL(clicked()), SLOT(on_click_dot()));
        connect(imagine, SIGNAL(clicked()), SLOT(on_click_imagine()));
        connect(add, SIGNAL(clicked()), SLOT(on_click_add()));
        connect(subtract, SIGNAL(clicked()), SLOT(on_click_subtract()));
        connect(multiply, SIGNAL(clicked()), SLOT(on_click_multiply()));
        connect(divide, SIGNAL(clicked()), SLOT(on_click_divide()));
        connect(equal, SIGNAL(clicked()), SLOT(on_click_equal()));
        connect(del, SIGNAL(clicked()), SLOT(on_click_del()));
        connect(clear, SIGNAL(clicked()), SLOT(on_click_clear()));
        connect(clear_entry, SIGNAL(clicked()), SLOT(on_click_clear_entry()));
        connect(left_bracket, SIGNAL(clicked()), SLOT(on_click_left_bracket()));
        connect(right_bracket, SIGNAL(clicked()), SLOT(on_click_right_bracket()));
        connect(sqroot, SIGNAL(clicked()), SLOT(on_click_sqrt()));
        connect(square, SIGNAL(clicked()), SLOT(on_click_square()));
        connect(percent, SIGNAL(clicked()), SLOT(on_click_percent()));
    }
    //Finished
    show();
}

void Calculator::keyPressEvent(QKeyEvent* event) {
    switch(event->key()) {
        case Qt::Key_0:
            on_click_num0();
            break;
        case Qt::Key_1:
            on_click_num1();
            break;
        case Qt::Key_2:
            on_click_num2();
            break;
        case Qt::Key_3:
            on_click_num3();
            break;
        case Qt::Key_4:
            on_click_num4();
            break;
        case Qt::Key_5:
            on_click_num5();
            break;
        case Qt::Key_6:
            on_click_num6();
            break;
        case Qt::Key_7:
            on_click_num7();
            break;
        case Qt::Key_8:
            on_click_num8();
            break;
        case Qt::Key_9:
            on_click_num9();
            break;
        case Qt::Key_Period:
            on_click_dot();
            break;
        case Qt::Key_I:
            on_click_imagine();
            break;
        case Qt::Key_Plus:
            on_click_add();
            break;
        case Qt::Key_Minus:
            on_click_subtract();
            break;
        case Qt::Key_Asterisk:
            on_click_multiply();
            break;
        case Qt::Key_Slash:
            on_click_divide();
            break;
        case Qt::Key_Equal:
        case Qt::Key_Return:
            on_click_equal();
            break;
        case Qt::Key_Backspace:
            on_click_del();
            break;
        case Qt::Key_ParenLeft:
            on_click_left_bracket();
            break;
        case Qt::Key_ParenRight:
            on_click_right_bracket();
            break;
        case Qt::Key_Percent:
            on_click_percent();
            break;
    }
}
////////////////////////////////////// On click buttons //////////////////////////////
void Calculator::on_click_num0() {
    editor_bottom->setText(editor_bottom->text() + num0->text());
}

void Calculator::on_click_num1() {
    editor_bottom->setText(editor_bottom->text() + num1->text());
}

void Calculator::on_click_num2() {
    editor_bottom->setText(editor_bottom->text() + num2->text());
}

void Calculator::on_click_num3() {
    editor_bottom->setText(editor_bottom->text() + num3->text());
}

void Calculator::on_click_num4() {
    editor_bottom->setText(editor_bottom->text() + num4->text());
}

void Calculator::on_click_num5() {
    editor_bottom->setText(editor_bottom->text() + num5->text());
}

void Calculator::on_click_num6() {
    editor_bottom->setText(editor_bottom->text() + num6->text());
}

void Calculator::on_click_num7() {
    editor_bottom->setText(editor_bottom->text() + num7->text());
}

void Calculator::on_click_num8() {
    editor_bottom->setText(editor_bottom->text() + num8->text());
}

void Calculator::on_click_num9() {
    editor_bottom->setText(editor_bottom->text() + num9->text());
}

void Calculator::on_click_dot() {
    editor_bottom->setText(editor_bottom->text() + dot->text());
}

void Calculator::on_click_imagine() {
    editor_bottom->setText(editor_bottom->text() + imagine->text());
}

void Calculator::on_click_add() {
    editor_bottom->setText(editor_bottom->text() + add->text());
}

void Calculator::on_click_subtract() {
    editor_bottom->setText(editor_bottom->text() + subtract->text());
}

void Calculator::on_click_multiply() {
    auto str = editor_bottom->text().toStdWString();
    str += L'×';
    editor_bottom->setText(QString::fromStdWString(str));
}

void Calculator::on_click_divide() {
    editor_bottom->setText(editor_bottom->text() + divide->text());
}

void Calculator::on_click_equal() {
    ExprReader reader(editor_bottom->text().toStdWString());
    auto result = reader.calc();
    editor_top->setText(editor_bottom->text());
    editor_bottom->setText(QString::fromStdString(result->toString()));
}

void Calculator::on_click_del() {
    int size = editor_bottom->text().size();
    if(size != 0)
        editor_bottom->setText(editor_bottom->text().left(size - 1));
}

void Calculator::on_click_clear() {
    editor_top->setText("");
    editor_bottom->setText("");
}

void Calculator::on_click_clear_entry() {
    editor_bottom->setText("");
}

void Calculator::on_click_left_bracket() {
    editor_bottom->setText(editor_bottom->text() + left_bracket->text());
}

void Calculator::on_click_right_bracket() {
    editor_bottom->setText(editor_bottom->text() + right_bracket->text());
}

void Calculator::on_click_sqrt() {
    ExprReader reader(editor_bottom->text().toStdWString());
    auto result = reader.calc();
    auto temp = *result * *result;
    editor_top->setText(editor_bottom->text());
    if(temp != nullptr) {
        *result << *temp;
        editor_bottom->setText(QString::fromStdString(result->toString()));
    }
    delete result;
}

void Calculator::on_click_square() {
    ExprReader reader(editor_bottom->text().toStdWString());
    auto result = reader.calc();
    auto temp = sqrt(*result);
    editor_top->setText(editor_bottom->text());
    if(temp != nullptr) {
        *result << *temp;
        editor_bottom->setText(QString::fromStdString(result->toString()));
    }
    delete result;
}

void Calculator::on_click_percent() {
    editor_bottom->setText(editor_bottom->text() + percent->text());
}