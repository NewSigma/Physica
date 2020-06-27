/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#include "Physica/Gui/LongtonAnt.h"

/*!
 * Inspired by interview question 16.22 on leetcode.
 * This is a kind of cellular automata.
 *
 * Optimize: QScatterSeries is implemented using QVector, use QList may improve performance.
 */
void LongtonAnt::handle() {
    int index;
    for(index = 0; index < series0->count(); ++index)
        if(series0->at(index).x() == r && series0->at(index).y() == c)
            break;

    if (index == series0->count()) {
        type = (type + 3) % 4;
        series0->append(r, c);
    } else {
        type = (type + 5) % 4;
        series0->remove(index);
    }

    switch(type) {
        case 0:
            ++r;
            break;
        case 1:
            ++c;
            break;
        case 2:
            --r;
            break;
        case 3:
            --c;
            break;
        default:;
    }

    X = -X < r - 1 ? X : 1 - r;
    X = X > r + 1 ? X : r + 1;
    Y = -Y < c - 1 ? Y : 1 - c;
    Y = Y > c + 1 ? Y : c + 1;
}

LongtonAnt::LongtonAnt() : QChartView(new QChart()), r(0), c(0), type(1), X(1), Y(1) {
    setAttribute(Qt::WA_DeleteOnClose);
    resize(650,650);

    series0 = new QScatterSeries(this);
    series0->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
    series0->setMarkerSize(10.0);
    startTimer(50);

    setRenderHint(QPainter::Antialiasing);
    chart()->addSeries(series0);

    chart()->setTitle("Longton's ant");
    chart()->createDefaultAxes();
    chart()->setDropShadowEnabled(false);

    show();
}

void LongtonAnt::timerEvent(QTimerEvent* event) {
    Q_UNUSED(event)
    handle();
    chart()->axes()[0]->setRange(QVariant(-X), QVariant(X));
    chart()->axes()[1]->setRange(QVariant(-Y), QVariant(Y));
}