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
#ifndef PHYSICA_LONGTONANT_H
#define PHYSICA_LONGTONANT_H

#include <set>
#include <QtCharts/QtCharts>

using namespace std;

namespace Physica::Gui {
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
}

#endif
