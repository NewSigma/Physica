/*
 * Copyright 2021 WeiBo He.
 *
 * This file is part of Physica.
 *
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
#include <iostream>
#include <QtWidgets/QApplication>
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Gui;
using namespace Physica::Utils;

int main(int argc, char** argv) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution phi{};
    Array<double> arr(5000);
    for (int i = 0; i < 5000; ++i)
        arr[i] = phi(gen);
    
    QApplication app(argc, argv);
    Plot* plot = new Plot();
    plot->hist(arr, 100);   
    plot->show();
    return QApplication::exec();
}
