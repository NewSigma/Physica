#include "Physica/Core/Math/Calculus/Differential.h"

namespace Physica::Core {
    Differential::Differential(TreeFunction func, MultiScalar at, MultiScalar stepSize)
            : func(std::move(func)), at(std::move(at)), stepSize(std::move(stepSize)) {}
    /*!
     * Optimize: if \at is much larger than \stepsize, the result will be 0. May be use talor series
     * and expend the function to the first order.
     */
    MultiScalar Differential::solve(DifferentialMethod method) {
        MultiScalar result;
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