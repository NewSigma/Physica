#ifndef PHYSICA_C_CALCULATOR_H
#define PHYSICA_C_CALCULATOR_H

#include <QtWidgets/QLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLineEdit>

namespace Physica::Gui {
    class Calculator : public QWidget {
    Q_OBJECT
    private:
        QVBoxLayout* default_layout;
        QLineEdit* editor_top;
        QLineEdit* editor_bottom;

        QGridLayout* layout;
        QPushButton* num0;
        QPushButton* num1;
        QPushButton* num2;
        QPushButton* num3;
        QPushButton* num4;
        QPushButton* num5;
        QPushButton* num6;
        QPushButton* num7;
        QPushButton* num8;
        QPushButton* num9;
        QPushButton* dot;
        QPushButton* imagine;
        QPushButton* add;
        QPushButton* subtract;
        QPushButton* multiply;
        QPushButton* divide;
        QPushButton* equal;
        QPushButton* del;
        QPushButton* clear;
        QPushButton* clear_entry;
        QPushButton* left_bracket;
        QPushButton* right_bracket;
        QPushButton* sqroot;
        QPushButton* square;
        QPushButton* percent;
    public:
        Calculator();
    private:
        void keyPressEvent(QKeyEvent* event) override;
    private Q_SLOTS:
        void on_click_num0();
        void on_click_num1();
        void on_click_num2();
        void on_click_num3();
        void on_click_num4();
        void on_click_num5();
        void on_click_num6();
        void on_click_num7();
        void on_click_num8();
        void on_click_num9();
        void on_click_dot();
        void on_click_imagine();
        void on_click_add();
        void on_click_subtract();
        void on_click_multiply();
        void on_click_divide();
        void on_click_equal();
        void on_click_del();
        void on_click_clear();
        void on_click_clear_entry();
        void on_click_left_bracket();
        void on_click_right_bracket();
        void on_click_sqrt();
        void on_click_square();
        void on_click_percent();
    };
}

#endif