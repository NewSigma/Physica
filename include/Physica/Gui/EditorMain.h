/*
 * Copyright 2020 WeiBo He.
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
