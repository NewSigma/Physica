/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtGui/QPainter>
#include <QtGui/QTextBlock>
#include "Gui/Header/LineNumberArea.h"
#include "Gui/Header/EditorMain.h"
#include <QScrollBar>

LineNumberArea::LineNumberArea(EditorMain* parent) : QWidget(parent), editor(parent) {
    setFont(editor->defaultFont);
}

void LineNumberArea::paintEvent(QPaintEvent* event) {
    QPainter painter(this);
    painter.fillRect(rect(), Qt::lightGray);

    auto rootFrame = editor->document()->rootFrame();
    auto frameIterator = rootFrame->begin();

    //Translate the coordinate system parallel with verticalScrollBar. And the y is the y-coordinate of old origin.
    double y = - editor->verticalScrollBar()->value();
    int index = 0;
    for(; frameIterator != rootFrame->end() && frameIterator.currentFrame() == nullptr; ++frameIterator) {
        qreal blockHeight = frameIterator.currentBlock().layout()->boundingRect().height();
        ++index;
        if(y + blockHeight > 0)
            break;
        y += blockHeight;
    }

    painter.setPen(Qt::black);
    for(; frameIterator != rootFrame->end() && frameIterator.currentFrame() == nullptr; ++frameIterator) {
        auto currentBlockRect = frameIterator.currentBlock().layout()->boundingRect();
        painter.drawText(0, static_cast<int>(y + (currentBlockRect.height() - fontMetrics().height()) / 2) + 5
                , event->rect().width(), fontMetrics().height()
                , Qt::AlignLeft, QString::number(index));
        y += currentBlockRect.height();
        if(y > editor->height())
            break;
        ++index;
    }
}