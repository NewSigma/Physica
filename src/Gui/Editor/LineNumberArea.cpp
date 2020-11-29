/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <QtGui/QPainter>
#include <QtGui/QTextBlock>
#include <QScrollBar>
#include "Physica/Gui/Editor/LineNumberArea.h"
#include "Physica/Gui/Editor/EditorMain.h"

namespace Physica::Gui {
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
}
