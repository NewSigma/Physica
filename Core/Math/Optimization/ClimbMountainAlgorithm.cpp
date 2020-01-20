#include "../../Header/ClimbMountainAlgorithm.h"
#include <iostream>
/*
 * The climb mountain algorithm should only be used in simple situation external point to get a global best solution.
 *
 * stepSize: Original step size that depend solving time we will use and whether can we get a global solution or not.
 * minStep: Final step size that depend the precision of our result.
 *
 * Operators that need memory free : func(arg)
 *
 * Usage: new ClimbMountainAlgorithm(args)->getExtremal();
 *
 * Copyright (c) 2019 NewSigma@163.com.All rights reserved.
 */
//TODO Optional: Auto step-size
ClimbMountainAlgorithm::ClimbMountainAlgorithm(RealNumber* f(RealNumber*), RealNumber* init, RealNumber* size) {
	func = f;
	x_initial = init;
	stepSize = size;
	minStep = new RealNumber(new unsigned char[1]{1}, 1, -16, true);
}

ClimbMountainAlgorithm::~ClimbMountainAlgorithm() {
    delete minStep;
}

void ClimbMountainAlgorithm::getExtremal() {
	RealNumber* y_initial = func(x_initial);
    auto x = new RealNumber();
    auto y = new RealNumber();
    auto x_last = new RealNumber();
    auto y_last = new RealNumber();
    auto x_result = new RealNumber();
    auto y_result = new RealNumber();
	bool positiveUsable = false;
	bool negativeUsable = false;

    *x_last = *x_initial;
    *x_last += *stepSize;
    *y_last = *func(x_last);
    if (*y_last > *y_initial) {
        positiveUsable = true;
        while (true) {
            *x = *x_last;
            *x += *stepSize;
            *y = *func(x);
            if (*y > *y_last) {
                *x_last = *x;
                *y_last = *y;
            }
            else {
                stepSize->power -= 1;
                if (*stepSize < *minStep) {
                    *x_result = *x_last;
                    *y_result = *y_last;
                    break;
                }
            }
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
                stepSize->power -= 1;
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
    std::cout << "[ClimbMountainAlgorithm] ";
    if (positiveUsable || negativeUsable) {
        std::cout << "External got: (";
        x_result->print();
        std::cout << ", ";
        y_result->print();
        std::cout << ")" <<  std::endl;
    }
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

RealNumber* ClimbMountainAlgorithm::getMinStep() {
	return minStep;
}

void ClimbMountainAlgorithm::setMinStep(RealNumber* step) {
    delete minStep;
	minStep = step;
}