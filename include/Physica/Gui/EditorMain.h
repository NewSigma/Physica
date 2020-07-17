#ifndef PHYSICA_EDITORMAIN_H
#define PHYSICA_EDITORMAIN_H

#include "QTextEdit"

namespace Physica::Gui {
    class LineNumberArea;

    class EditorMain : public QTextEdit {
        LineNumberArea* lineNumberArea;
        //Should be put into settings.
        QFont defaultFont;
    public:
        explicit EditorMain(QWidget* parent);
    protected:
        void resizeEvent(QResizeEvent* event) override;
    private:
        int lineNumberAreaWidth() const;
        void doHighLight();
    private slots:
        void updateLineNumberAreaWidth();
        void onCursorPositionChanged();

        friend class LineNumberArea;
    };
}

#endif
