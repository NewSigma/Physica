#ifndef PHYSICA_EDITORMAIN_H
#define PHYSICA_EDITORMAIN_H

#include "QTextEdit"

class EditorMain : public QTextEdit {
public:
    explicit EditorMain(QWidget* parent);
private:
    int highLight;
    int lineNumberAreaWidth();
    void doHighLight();
    //void paintLineNumber(QPaintEvent *event);
private slots:
    void updateLineNumberAreaWidth();
    void onCursorPositionChanged();
};


#endif
