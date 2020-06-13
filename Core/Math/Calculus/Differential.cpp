#include <Core/Header/Differential.h>
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    Differential::Differential(Function func, Numerical at, Numerical stepSize)
            : func(std::move(func)), at(std::move(at)), stepSize(std::move(stepSize)) {}
    /*!
     * Optimize: if \at is much larger than \stepsize, the result will be 0. May be use talor series
     * and expend the function to the first order.
     */
    Numerical Differential::solve(DifferentialMethod method) {
        Numerical result;
        switch(method) {
            case DoublePoint:
                result = (func(at + stepSize) - func(at - stepSize)) / (stepSize << 1);
                break;
            case Forward:
                result = (func(at + stepSize) - func(at)) / stepSize;
                break;
            case Backward:
                result = (func(at) - func(at - stepSize)) / stepSize;
                break;
        }
        return result;
    }
}