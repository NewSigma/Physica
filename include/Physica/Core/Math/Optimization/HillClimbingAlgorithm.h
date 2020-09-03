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
#ifndef PHYSICA_HILLCLIMBINGALGORITHM_H
#define PHYSICA_HILLCLIMBINGALGORITHM_H

#include "Physica/Core/Math/Calculus/Function/VectorFunction/VectorFunction.h"
#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica::Core {
    /*!
     * @class HillClimbingAlgorithm tries to get the maximum of a function.
     * If you want the minimum of a function f(x), simply calculate the maximum of -f(x).
     */
    template<size_t dim, ScalarType type>
    class HillClimbingAlgorithm {
        static_assert(dim > 0, "dim must be at least 1.");
    };

    template<ScalarType type>
    class HillClimbingAlgorithm<1, type> {
    public:
        enum State {
            Unavailable,
            Ready,
            OverFlow
        };
    private:
        const VectorFunction<type, false> func;
        Scalar<type, false> x_initial;
        /*!
         * \minStep is the minimum step size that depends the precision of our result.
         * Set a positive value to search from the positive axe of x,
         * or a negative value to search from negative value.
         */
        Scalar<type, false> minStep;
        State state;
    public:
        HillClimbingAlgorithm(const VectorFunction<type, false>& func
                              , const Scalar<type, false> initial
                              , const Scalar<type, false> minStep);
        HillClimbingAlgorithm(const HillClimbingAlgorithm& alg);
        HillClimbingAlgorithm(HillClimbingAlgorithm& alg);
        ~HillClimbingAlgorithm() = default;
        /* Operations */
        Point<2, type, false> solve() const;
        /* Getters */
        [[nodiscard]] Scalar<type, false> getMinStep() const { return minStep; }
        [[nodiscard]] State getState() const { return state; }
        /* Setters */
        void setMinStep(const Scalar<type, false>& s) { minStep = s; }
    };
}

#include "HillClimbingAlgorithmImpl.h"

#endif