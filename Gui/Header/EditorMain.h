#ifndef PHYSICA_EDITORMAIN_H
#define PHYSICA_EDITORMAIN_H

#include "QTextEdit"

class EditorMain : public QTextEdit {
public:
    EditorMain(QWidget* parent);
private:
    int highLight;
    int lineNumberAreaWidth();
    void paintLineNumber(QPaintEvent *event);
private slots:
    void updateLineNumberAreaWidth();
    void doHighLight();
};


#endif
