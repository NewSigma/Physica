#include <Physica/Core/Integrate.h>
#include "Physica/Core/Scalar.h"

namespace Physica::Core {
    Integrate::Integrate(Function func, Scalar from, Scalar to, Scalar stepSize)
            : func(std::move(func)), from(std::move(from)), to(std::move(to)), stepSize(std::move(stepSize)) {}
    /*!
     * Optimize: if \at is much larger than \stepsize, the result will be 0. May be use talor series
     * and expend the function to the first order.
     */
    Scalar Integrate::solve(IntegrateMethod method) {
        Scalar result(BasicConst::getInstance().get_0());
        switch(method) {
            case Rectangular: {
                Scalar start(from);
                while(start < to) {
                    result += func(start);
                    start += stepSize;
                }
                result *= stepSize;
            }
                break;
            case Ladder: {
                result += ((func(from) + func(to)) >> 1);
                Scalar start(from + stepSize);
                while(start < to) {
                    result += func(start);
                    start += stepSize;
                }
                result *= stepSize;
            }
                break;
            case Simpson: {
                result += func(from) + func(to);
                Scalar odd(BasicConst::getInstance().get_0());
                Scalar even(BasicConst::getInstance().get_0());
                bool b = true;
                Scalar start = from + stepSize;
                while(start < to) {
                    Scalar& toChange = b ? odd : even;
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