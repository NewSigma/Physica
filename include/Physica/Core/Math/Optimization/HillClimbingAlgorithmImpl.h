/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

namespace Physica::Core {
    template<ScalarOption option>
    HillClimbingAlgorithm<1, option>::HillClimbingAlgorithm(const VectorFunction<option, false>& func
            , const ScalarType initial, const ScalarType minStep)
            : func(func), x_initial(initial), minStep(minStep) {
        state = minStep.isZero() ? Unavailable : Ready;
    }

    template<ScalarOption option>
    HillClimbingAlgorithm<1, option>::HillClimbingAlgorithm(const HillClimbingAlgorithm<1, option>& alg)
            : func(alg.func), x_initial(alg.x_initial), minStep(alg.minStep), state(alg.state) {}

    template<ScalarOption option>
    HillClimbingAlgorithm<1, option>::HillClimbingAlgorithm(HillClimbingAlgorithm<1, option>& alg)
            : func(std::move(alg.func))
            , x_initial(std::move(alg.x_initial))
            , minStep(std::move(alg.minStep))
            , state(alg.state) {}

    template<ScalarOption option>
    Point<2, typename HillClimbingAlgorithm<1, option>::ScalarType> HillClimbingAlgorithm<1, option>::solve() const {
        if(state != Ready)
            return func(x_initial);

        const ScalarType y_initial = func(x_initial);
        ScalarType x_last(x_initial);
        ScalarType y_last(y_initial);
        ScalarType x = x_initial + minStep;
        ScalarType y = func(x_last);
        ScalarType x_result{};
        ScalarType y_result{};
        ScalarType stepSize(minStep);
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
        return Point<2, ScalarType>{x_result, y_result};
    }
}
