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
#ifndef PHYSICA_HILLCLIMBINGALGORITHMIMPL_H
#define PHYSICA_HILLCLIMBINGALGORITHMIMPL_H

namespace Physica::Core {
    template<ScalarType type>
    HillClimbingAlgorithm<1, type>::HillClimbingAlgorithm(const VectorFunction<type, false>& func
            , const Scalar<type, false> initial, const Scalar<type, false> minStep)
            : func(func), x_initial(initial), minStep(minStep) {
        state = minStep.isZero() ? Unavailable : Ready;
    }

    template<ScalarType type>
    HillClimbingAlgorithm<1, type>::HillClimbingAlgorithm(const HillClimbingAlgorithm<1, type>& alg)
            : func(alg.func), x_initial(alg.x_initial), minStep(alg.minStep), state(alg.state) {}

    template<ScalarType type>
    HillClimbingAlgorithm<1, type>::HillClimbingAlgorithm(HillClimbingAlgorithm<1, type>& alg)
            : func(std::move(alg.func))
            , x_initial(std::move(alg.x_initial))
            , minStep(std::move(alg.minStep))
            , state(alg.state) {}

    template<ScalarType type>
    Point<2, type, false> HillClimbingAlgorithm<1, type>::solve() const {
        if(state != Ready)
            return func(x_initial);

        const Scalar<type, false> y_initial = func(x_initial);
        Scalar<type, false> x_last(x_initial);
        Scalar<type, false> y_last(y_initial);
        Scalar<type, false> x = x_initial + minStep;
        Scalar<type, false> y = func(x_last);
        Scalar<type, false> x_result{};
        Scalar<type, false> y_result{};
        Scalar<type, false> stepSize(minStep);
        /*
         * We use changeable stepSize, if y > y_last, stepMultiple += 1, otherwise stepMultiple -= 1.
         * The current step size is minStep * (stepMultiple + 1).
         */
        unsigned char stepMultiple = 0;
        //If the minStep is too long, some extremal points may be ignored.
        while (true) {
            if (y > y_last) {
                x_last = x;
                y_last = y;
                stepSize <<= 1;
                ++stepMultiple;
                if(stepMultiple == UCHAR_MAX) {
                    state = OverFlow;
                    return y_last;
                }
            }
            else {
                stepSize >>= 1;
                --stepMultiple;
                if(stepMultiple == 0) {
                    x_result = x_last;
                    y_result = y_last;
                    break;
                }
            }
            x = x_last + stepSize;
            y = func(x);
        }
        return Point<2, type, false>(x_result, y_result);
    }
}

#endif