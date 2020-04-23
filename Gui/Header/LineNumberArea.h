/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LINENUMBERAREA_H
#define PHYSICA_LINENUMBERAREA_H

#include "QWidget"

class EditorMain;

class LineNumberArea : public QWidget {
    const EditorMain* const editor;
public:
    explicit LineNumberArea(EditorMain* parent);
protected:
    void paintEvent(QPaintEvent* event) override;

    Q_DISABLE_COPY_MOVE(LineNumberArea)
};

#endif