/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtGui/QPainter>
#include <QtGui/QTextBlock>
#include <QPaintEvent>
#include <QScrollBar>
#include "Physica/Gui/Plot.h"
#include "Physica/Gui/EditorMain.h"
#include "Physica/Gui/LineNumberArea.h"

namespace Physica::Gui {
    EditorMain::EditorMain(QWidget* parent) : QTextEdit(parent), defaultFont(QFont("DejaVu Sans Mono", 14)) {
        lineNumberArea = new LineNumberArea(this);

        document()->setDefaultFont(defaultFont);

        updateLineNumberAreaWidth();
        doHighLight();

        connect(document(), &QTextDocument::blockCountChanged, this, &EditorMain::updateLineNumberAreaWidth);
        connect(this, &EditorMain::cursorPositionChanged, this, &EditorMain::onCursorPositionChanged);
        connect(this, SIGNAL(cursorPositionChanged()), lineNumberArea, SLOT(update()));
        connect(verticalScrollBar(), SIGNAL(valueChanged(int)), lineNumberArea, SLOT(update()));
    }

    void EditorMain::resizeEvent(QResizeEvent* event) {
        lineNumberArea->setFixedHeight(height());
        QTextEdit::resizeEvent(event);
    }

    int EditorMain::lineNumberAreaWidth() const {
        int digits = 1;
        int max = qMax(1, document()->blockCount());
        while (max >= 10) {
            max /= 10;
            ++digits;
        }
        return 80 + fontMetrics().horizontalAdvance(QLatin1Char('9')) * digits;
    }

    void EditorMain::doHighLight() {
        QList<QTextEdit::ExtraSelection> extraSelections;
        if (!isReadOnly()) {
            QTextEdit::ExtraSelection selection{};

            selection.format.setBackground(QColor(0xDCDCDC));
            selection.format.setProperty(QTextFormat::FullWidthSelection, true);
            selection.cursor = textCursor();
            selection.cursor.clearSelection();
            extraSelections.append(selection);
        }
        setExtraSelections(extraSelections);
    }

    void EditorMain::updateLineNumberAreaWidth() {
        int w = lineNumberAreaWidth();
        setViewportMargins(w, 0, 0, 0);
        lineNumberArea->setFixedWidth(w);
    }

    void EditorMain::onCursorPositionChanged() {
        doHighLight();
    }
}
