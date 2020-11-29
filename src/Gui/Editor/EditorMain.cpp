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
#include <QPaintEvent>
#include <QScrollBar>
#include "Physica/Gui/Plot/Plot.h"
#include "Physica/Gui/Editor/EditorMain.h"
#include "Physica/Gui/Editor/LineNumberArea.h"

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
