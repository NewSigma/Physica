#include <Core/Header/Integrate.h>
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    Integrate::Integrate(Function func, Numerical from, Numerical to, Numerical stepSize)
            : func(std::move(func)), from(std::move(from)), to(std::move(to)), stepSize(std::move(stepSize)) {}
    /*!
     * Optimize: if \at is much larger than \stepsize, the result will be 0. May be use talor series
     * and expend the function to the first order.
     */
    Numerical Integrate::solve(IntegrateMethod method) {
        Numerical result(BasicConst::getInstance().get_0());
        switch(method) {
            case Rectangular: {
                Numerical start(from);
                while(start < to) {
                    result += func(start);
                    start += stepSize;
                }
                result *= stepSize;
            }
                break;
            case Ladder: {
                result += ((func(from) + func(to)) >> 1);
                Numerical start(from + stepSize);
                while(start < to) {
                    result += func(start);
                    start += stepSize;
                }
                result *= stepSize;
            }
                break;
            case Simpson: {
                result += func(from) + func(to);
                Numerical odd(BasicConst::getInstance().get_0());
                Numerical even(BasicConst::getInstance().get_0());
                bool b = true;
                Numerical start = from + stepSize;
                while(start < to) {
                    Numerical& toChange = b ? odd : even;
                    b = !b;
                    toChange += func(start);
                    start += stepSize;
                }
                odd <<= 2;
                even <<= 1;
                result += odd + even;
                result *= stepSize;
                result /= BasicConst::getInstance().get_3();
            }
                break;
            case Simpson_3_8:
                break;
            case Bode:
                break;
        }
        return result;
    }
}