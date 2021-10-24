/*
 * Copyright 2021 WeiBo He.
 *
 * This file is part of Physica.
 *
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

#include <cassert>
#include <functional>
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    template<class ScalarType> class ScalarBase;
    template<class VectorType> class RValueVector;
}

namespace Physica::Core::Math {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.121
     */
    template<class ScalarType, class Function, class VectorType>
    class ConjugateGradient {
        constexpr static double c1 = 1E-4;
        constexpr static double c2 = 0.1;
        constexpr static size_t maxIteration = 100000;

        Function func;
        VectorType x;
        ScalarType epsilon;
        ScalarType maxStepSize;
        ScalarType diffStep;
    public:
        ConjugateGradient(Function func_,
                          const RValueVector<VectorType>& x_,
                          const ScalarBase<ScalarType>& epsilon_,
                          const ScalarBase<ScalarType>& maxStepSize_ = ScalarType(1),
                          const ScalarBase<ScalarType>& diffStep_ = ScalarType(1E-6));
        ~ConjugateGradient() = default;
        /* Operations */
        ScalarType compute();
        /* Getters */
        [[nodiscard]] const VectorType& getX() const noexcept { return x; }
    private:
        [[nodiscard]] VectorType getGradient(const VectorType& at) const;
        [[nodiscard]] ScalarType lineSearch(const VectorType& gradient, const VectorType& direction);
        [[nodiscard]] ScalarType zoom(const VectorType& gradient, const VectorType& direction, ScalarType step1, ScalarType step2);
        [[nodiscard]] ScalarType interpolate(const VectorType& direction, const ScalarType& step1, const ScalarType& step2);
    };

    template<class ScalarType, class Function, class VectorType>
    ConjugateGradient<ScalarType, Function, VectorType>::ConjugateGradient(Function func_,
                                                                           const RValueVector<VectorType>& x_,
                                                                           const ScalarBase<ScalarType>& epsilon_,
                                                                           const ScalarBase<ScalarType>& maxStepSize_,
                                                                           const ScalarBase<ScalarType>& diffStep_)
            : func(func_)
            , x(x_.getDerived())
            , epsilon(epsilon_.getDerived())
            , maxStepSize(maxStepSize_.getDerived())
            , diffStep(diffStep_.getDerived()) {
        assert(maxStepSize.isPositive());
        assert(diffStep.isPositive());
    }

    template<class ScalarType, class Function, class VectorType>
    ScalarType ConjugateGradient<ScalarType, Function, VectorType>::compute() {
        VectorType g = getGradient(x);
        VectorType direction = -g;
        ScalarType error = g.squaredNorm();

        size_t ite = 0;
        const size_t N = x.getLength();
        while (error > epsilon) {
            const ScalarType factor = lineSearch(g, direction);
            x += factor * direction;

            const VectorType g1 = getGradient(x);
            const ScalarType new_error = g1.squaredNorm();
            if (++ite == N) {
                direction = -g1;
                ite = 0;
            }
            else {
                const ScalarType alpha = (new_error - g * g1) / error;
                direction = alpha * direction - g1;
            }
            g = g1;
            error = new_error;
        }
        return func(std::cref(x));
    }

    template<class ScalarType, class Function, class VectorType>
    VectorType ConjugateGradient<ScalarType, Function, VectorType>::getGradient(const VectorType& at) const {
        const size_t length = at.getLength();
        VectorType result(length);
        VectorType copy = at;
        for (size_t i = 0; i < length; ++i) {
            const ScalarType buffer = copy[i];
            result[i] = Differential<ScalarType>::doublePoint([&](ScalarType alpha) {
                copy[i] = alpha;
                return func(std::cref(copy));
            }, buffer, diffStep);
            copy[i] = buffer;
        }
        return result;
    }

    template<class ScalarType, class Function, class VectorType>
    ScalarType ConjugateGradient<ScalarType, Function, VectorType>::lineSearch(
            const VectorType& gradient,
            const VectorType& direction) {
        ScalarType step_lower = ScalarType::Zero();
        const ScalarType step_upper = maxStepSize;
        ScalarType step = maxStepSize / ScalarType::Two();

        const ScalarType factor1 = ScalarType(c1);
        const ScalarType factor2 = ScalarType(c2);
        const ScalarType phi_0 = func(std::cref(x));
        const ScalarType diff_phi_0 = gradient * direction;
        assert(diff_phi_0.isNegative());

        ScalarType last_y = phi_0;
        VectorType x1 = VectorType(x.getLength());
        size_t i = 0;
        while (true) {
            x1 = x + step * direction;
            const ScalarType y = func(std::cref(x1));
            const bool violatesWolfe = y > (phi_0 + factor1 * step * diff_phi_0);
            const bool isIncreased = (y >= last_y) && (i > 0);
            if (violatesWolfe || isIncreased)
                return zoom(gradient, direction, step_lower, step);
            
            const ScalarType diff_phi = getGradient(x1) * direction;
            if (abs(diff_phi) <= -factor2 * diff_phi_0)
                return step;
            if (!diff_phi.isNegative())
                return zoom(gradient, direction, step, step_lower);
            
            step_lower = step;
            step = (step_lower + step_upper) / ScalarType::Two();
            if (abs(step_upper - step) < std::numeric_limits<ScalarType>::epsilon())
                return step_upper;
            if (++i > maxIteration)
                throw BadConvergenceException();
        }
    }
    /**
     * \param step1 does not necessarily less than \param step2
     */
    template<class ScalarType, class Function, class VectorType>
    ScalarType ConjugateGradient<ScalarType, Function, VectorType>::zoom(
            const VectorType& gradient,
            const VectorType& direction,
            ScalarType step1,
            ScalarType step2) {
        const ScalarType phi_0 = func(std::cref(x));
        const ScalarType diff_phi_0 = gradient * direction;
        const ScalarType factor1 = ScalarType(c1);
        const ScalarType factor2 = ScalarType(c2);

        ScalarType last_y = func(VectorType(x + step1 * direction));
        VectorType x1 = VectorType(x.getLength());
        while (true) {
            const ScalarType step = interpolate(direction, step1, step2);

            x1 = x + step * direction;
            const ScalarType y = func(std::cref(x1));
            const bool violatesWolfe = y > (phi_0 + factor1 * step * diff_phi_0);
            const bool isIncreased = y >= last_y;
            if (violatesWolfe || isIncreased)
                step2 = step;
            else {
                const ScalarType diff_phi = getGradient(x1) * direction;
                if (abs(diff_phi) <= -factor2 * diff_phi_0)
                    return step;
                if (!(diff_phi * (step2 - step1)).isNegative())
                    step2 = step1;
                step1 = step;
                last_y = y;
            }

            if (abs(step2 - step1) < std::numeric_limits<ScalarType>::epsilon())
                return step1;
        }
    }

    template<class ScalarType, class Function, class VectorType>
    ScalarType ConjugateGradient<ScalarType, Function, VectorType>::interpolate(
            const VectorType& direction,
            const ScalarType& step1,
            const ScalarType& step2) {
        assert(step1 != step2);
        const VectorType x1 = x + step1 * direction;
        const VectorType x2 = x + step2 * direction;
        const ScalarType diff1 = getGradient(x1) * direction;
        const ScalarType diff2 = getGradient(x2) * direction;
        const ScalarType delta_step = step1 - step2;
        const ScalarType d1 = diff1 + diff2 - ScalarType(3) * (func(x1) - func(x2)) / delta_step;
        const ScalarType squared_d2 = square(d1) - diff1 * diff2;
        assert(!squared_d2.isNegative());
        ScalarType d2 = sqrt(squared_d2);
        if (delta_step.isNegative())
            d2.toOpposite();
        return step1 - delta_step * (diff1 + d2 - d1) / (diff1 - diff2 + ScalarType::Two() * d2);
    }
}
