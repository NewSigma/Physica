#ifndef PHYSICA_LONGTONANT_H
#define PHYSICA_LONGTONANT_H

#include <set>
#include "QtCharts"

using namespace std;

class LongtonAnt : public QChartView {
    Q_OBJECT
public:
    int r, c, X, Y;
    //D = 0, R = 1, U + 2, L = 3
    int type;
    QScatterSeries* series0;
    explicit LongtonAnt();
    void timerEvent(QTimerEvent* event) override;
    void handle();
};

#endif
