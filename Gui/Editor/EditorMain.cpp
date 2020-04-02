/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtGui/QPainter>
#include <QtGui/QTextBlock>
#include "../Header/EditorMain.h"
#include "QPaintEvent"

EditorMain::EditorMain(QWidget* parent) : QTextEdit(parent) {
    highLight = 0xDCDCDC;
    updateLineNumberAreaWidth();
    doHighLight();

    document()->rootFrame()->begin().currentBlock().setVisible(false);
    setReadOnly(true);

    connect(document(), &QTextDocument::blockCountChanged, this, &EditorMain::updateLineNumberAreaWidth);
    //connect(this, SIGNAL(updateRequest()), SLOT(updateLineNumberArea()));
    connect(this, &EditorMain::cursorPositionChanged, this, &EditorMain::onCursorPositionChanged);
}

int EditorMain::lineNumberAreaWidth() {
    int digits = 1;
    int max = qMax(1, document()->blockCount());
    while (max >= 10) {
        max /= 10;
        ++digits;
    }
    return 80 + fontMetrics().horizontalAdvance(QLatin1Char('9')) * digits;
}

void EditorMain::updateLineNumberAreaWidth() {
    setViewportMargins(lineNumberAreaWidth(), 0, 0, 0);
}

void EditorMain::doHighLight() {
    QList<QTextEdit::ExtraSelection> extraSelections;
    if (!isReadOnly()) {
        QTextEdit::ExtraSelection selection{};

        selection.format.setBackground(QColor(highLight));
        selection.format.setProperty(QTextFormat::FullWidthSelection, true);
        selection.cursor = textCursor();
        selection.cursor.clearSelection();
        extraSelections.append(selection);
    }
    setExtraSelections(extraSelections);
}

void EditorMain::onCursorPositionChanged() {
    setReadOnly(textCursor().currentFrame() == document()->rootFrame());
    doHighLight();
}
/*
void EditorMain::paintLineNumber(QPaintEvent *event) {
    QPainter painter(this);
    painter.fillRect(event->rect(), Qt::lightGray);

    QTextBlock block = document()->firstBlock();
    int blockNumber = block.blockNumber();
    int top = (int)blockBoundingGeometry(block).translated(contentOffset()).top();
    int bottom = top + (int)blockBoundingRect(block).height();

    while (block.isValid() && top <= event->rect().bottom()) {
        if (block.isVisible() && bottom >= event->rect().top()) {
            painter.setPen(Qt::black);
            painter.drawText(0, top, lineNumberAreaWidth(), fontMetrics().height(), Qt::AlignLeft, QString::number(blockNumber + 1));
        }

        block = block.next();
        top = bottom;
        bottom = top + (int)blockBoundingRect(block).height();
        ++blockNumber;
    }
}
*/