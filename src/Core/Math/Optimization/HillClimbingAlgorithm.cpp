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
#include <iostream>
#include "Physica/Core/Math/Optimization/HillClimbingAlgorithm.h"

namespace Physica::Core {

    //TODO Optional: Auto step-size
    HillClimbingAlgorithm::HillClimbingAlgorithm(MultiScalar* f(MultiScalar*), MultiScalar* init, MultiScalar* size) {
        func = f;
        x_initial = init;
        stepSize = size;
        minStep = new MultiScalar(BasicConst::getInstance().getStepSize());
    }

    HillClimbingAlgorithm::~HillClimbingAlgorithm() {
        delete minStep;
    }

    void HillClimbingAlgorithm::getExtremal() {
        MultiScalar* y_initial = func(x_initial);
        auto x_last = new MultiScalar(*x_initial);
        *x_last += *stepSize;
        auto y_last = func(x_last);
        auto x = new MultiScalar(*x_last);
        auto y = new MultiScalar(*stepSize);
        MultiScalar* x_result = nullptr;
        MultiScalar* y_result = nullptr;
        bool positiveUsable = false;
        bool negativeUsable = false;

        if (*y_last > *y_initial) {
            positiveUsable = true;
            while (true) {
                if (*y > *y_last) {
                    *x_last = *x;
                    *y_last = *y;
                }
                else {
                    *stepSize >> 1;
                    if (*stepSize < *minStep) {
                        x_result = x_last;
                        y_result = y_last;
                        break;
                    }
                }
                *x = *x_last;
                *x += *stepSize;
                *y = *func(x);
            }
        }

        *x_last = *x_initial;
        *x_last -= *stepSize;
        *y_last = *func(x_last);
        if (*y_last > *y_initial) {
            negativeUsable = true;
            while (true) {
                *x = *x_last;
                *x -= *stepSize;
                *y = *func(x);
                if (*y > *y_last) {
                    *x_last = *x;
                    *y_last = *y;
                }
                else {
                    *stepSize >> 1;
                    if (*stepSize < *minStep) {
                        if (positiveUsable && *y_last > *y_result) {
                            *x_result = *x_last;
                            *y_result = *y_last;
                        }
                        break;
                    }
                }
            }
        }
        //Print final result.
        std::cout << "[HillClimbingAlgorithm] ";
        if (positiveUsable || negativeUsable)
            std::cout << "External got: (" << *x_result << ", " << *y_result << ")" << std::endl;
        else
            std::cout << "Failed to fetch a external!";

        delete y_initial;
        delete x;
        delete y;
        delete x_last;
        delete y_last;
        delete x_result;
        delete y_result;
    }

    MultiScalar* HillClimbingAlgorithm::getMinStep() {
        return minStep;
    }

    void HillClimbingAlgorithm::setMinStep(MultiScalar* step) {
        delete minStep;
        minStep = step;
    }
}
