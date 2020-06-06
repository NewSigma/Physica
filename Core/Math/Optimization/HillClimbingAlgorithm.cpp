/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iostream>
#include "Core/Header/HillClimbingAlgorithm.h"
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    /*
     * The climb mountain algorithm should only be used in simple situation external point to get a global best solution.
     *
     * stepSize: Original step size that depend solving time we will use and whether can we get a global solution or not.
     * minStep: Final step size that depend the precision of our result.
     *
     * Operators that need memory free : func(arg)
     *
     * Usage: new HillClimbingAlgorithm(args)->getExtremal();
     */
    //TODO Optional: Auto step-size
    HillClimbingAlgorithm::HillClimbingAlgorithm(Numerical* f(Numerical*), Numerical* init, Numerical* size) {
        func = f;
        x_initial = init;
        stepSize = size;
        minStep = new Numerical(BasicConst::getInstance().getStepSize());
    }

    HillClimbingAlgorithm::~HillClimbingAlgorithm() {
        delete minStep;
    }

    void HillClimbingAlgorithm::getExtremal() {
        Numerical* y_initial = func(x_initial);
        auto x_last = new Numerical(*x_initial);
        *x_last += *stepSize;
        auto y_last = func(x_last);
        auto x = new Numerical(*x_last);
        auto y = new Numerical(*stepSize);
        Numerical* x_result = nullptr;
        Numerical* y_result = nullptr;
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

    Numerical* HillClimbingAlgorithm::getMinStep() {
        return minStep;
    }

    void HillClimbingAlgorithm::setMinStep(Numerical* step) {
        delete minStep;
        minStep = step;
    }
}
