/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.56-62
     */
    template<class ScalarType, size_t Dim>
    class LineSearch {
    public:
        using VectorType = Vector<ScalarType, Dim>;

        ScalarType maxStepSize;
        ScalarType decreaseCondNum;
        ScalarType curvatureCondNum;
        size_t maxIteration;
    public:
        LineSearch(ScalarType maxStepSize_, ScalarType decreaseCondNum_ = 1E-4, ScalarType curvatureCondNum_ = 0.1, size_t maxIteration_ = 100);
        LineSearch(const LineSearch&) = default;
        LineSearch(LineSearch&&) noexcept = default;
        ~LineSearch() = default;
        /* Operators */
        LineSearch& operator=(LineSearch obj) noexcept;
        /* Operations */
        template<class Functor, class GradFunctor>
        [[nodiscard]] ScalarType run(Functor func, GradFunctor grad, const VectorType& x, const VectorType& gradient, const VectorType& direction) const;
        void swap(LineSearch& __restrict obj) noexcept;
    private:
        template<class Functor, class GradFunctor>
        [[nodiscard]] ScalarType zoom(Functor func, GradFunctor grad, const VectorType& x, const VectorType& gradient, const VectorType& direction, ScalarType step1, ScalarType step2) const;
        template<class Functor, class GradFunctor>
        [[nodiscard]] ScalarType interpolate(Functor func, GradFunctor grad, const VectorType& x, const VectorType& direction, const ScalarType& step1, const ScalarType& step2) const;
    };

    template<class ScalarType, size_t Dim>
    LineSearch<ScalarType, Dim>::LineSearch(
                ScalarType maxStepSize_,
                ScalarType decreaseCondNum_,
                ScalarType curvatureCondNum_,
                size_t maxIteration_)
            : maxStepSize(maxStepSize_)
            , decreaseCondNum(decreaseCondNum_)
            , curvatureCondNum(curvatureCondNum_)
            , maxIteration(maxIteration_) {
        assert(maxStepSize.isPositive());
        assert(decreaseCondNum.isPositive());
        assert(decreaseCondNum < curvatureCondNum);
        assert(curvatureCondNum < ScalarType(1));
    }

    template<class ScalarType, size_t Dim>
    LineSearch<ScalarType, Dim>& LineSearch<ScalarType, Dim>::operator=(LineSearch<ScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    ScalarType LineSearch<ScalarType, Dim>::run(
            Functor func,
            GradFunctor grad,
            const VectorType& x,
            const VectorType& gradient,
            const VectorType& direction) const {
        if (gradient.squaredNorm() < std::numeric_limits<ScalarType>::min())
            return 0;
        ScalarType step_lower = ScalarType(0);
        const ScalarType step_upper = maxStepSize;
        ScalarType step = maxStepSize / ScalarType(2);

        const ScalarType phi_0 = func(x);
        const ScalarType diff_phi_0 = gradient * direction;
        assert(diff_phi_0.isNegative());

        ScalarType last_y = phi_0;
        VectorType x1 = VectorType(x.getLength());
        size_t i = 0;
        while (true) {
            x1 = x + step * direction;
            const ScalarType y = func(x1);
            const bool violatesWolfe = y > (phi_0 + decreaseCondNum * step * diff_phi_0);
            const bool isIncreased = (y >= last_y) && (i > 0);
            if (violatesWolfe || isIncreased)
                return zoom(func, grad, x, gradient, direction, step_lower, step);
            
            const ScalarType diff_phi = grad(x1) * direction;
            if (abs(diff_phi) <= -curvatureCondNum * diff_phi_0)
                return step;
            if (!diff_phi.isNegative())
                return zoom(func, grad, x, gradient, direction, step, step_lower);
            
            step_lower = step;
            step = (step_lower + step_upper) / ScalarType(2);
            if (abs(step_upper - step) < std::numeric_limits<ScalarType>::epsilon())
                return step_upper;
            if (++i > maxIteration)
                throw BadConvergenceException("Exceed max iteration of LineSearch");
        }
    }

    template<class ScalarType, size_t Dim>
    void LineSearch<ScalarType, Dim>::swap(LineSearch& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        maxStepSize.swap(obj.maxStepSize);
        decreaseCondNum.swap(obj.decreaseCondNum);
        curvatureCondNum.swap(obj.curvatureCondNum);
        std::swap(maxIteration, obj.maxIteration);
    }
    /**
     * \param step1 does not necessarily less than \param step2
     */
    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    ScalarType LineSearch<ScalarType, Dim>::zoom(
            Functor func,
            GradFunctor grad,
            const VectorType& x,
            const VectorType& gradient,
            const VectorType& direction,
            ScalarType step1,
            ScalarType step2) const {
        const ScalarType phi_0 = func(x);
        const ScalarType diff_phi_0 = gradient * direction;

        ScalarType last_y = func(VectorType(x + step1 * direction));
        VectorType x1 = VectorType(x.getLength());
        while (true) {
            const ScalarType step = interpolate(func, grad, x, direction, step1, step2);

            x1 = x + step * direction;
            const ScalarType y = func(x1);
            const bool violatesWolfe = y > (phi_0 + decreaseCondNum * step * diff_phi_0);
            const bool isIncreased = y >= last_y;
            if (violatesWolfe || isIncreased)
                step2 = step;
            else {
                const ScalarType diff_phi = grad(x1) * direction;
                if (abs(diff_phi) <= -curvatureCondNum * diff_phi_0)
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

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    ScalarType LineSearch<ScalarType, Dim>::interpolate(
            Functor func,
            GradFunctor grad,
            const VectorType& x,
            const VectorType& direction,
            const ScalarType& step1,
            const ScalarType& step2) const {
        assert(step1 != step2);
        const VectorType x1 = x + step1 * direction;
        const VectorType x2 = x + step2 * direction;
        const ScalarType diff1 = grad(x1) * direction;
        const ScalarType diff2 = grad(x2) * direction;
        const ScalarType delta_step = step1 - step2;
        const ScalarType d1 = diff1 + diff2 - ScalarType(3) * (func(x1) - func(x2)) / delta_step;
        const ScalarType squared_d2 = square(d1) - diff1 * diff2;
        assert(!squared_d2.isNegative());
        ScalarType d2 = sqrt(squared_d2);
        if (delta_step.isNegative())
            d2.toOpposite();
        return step1 - delta_step * (diff1 + d2 - d1) / (diff1 - diff2 + ScalarType(2) * d2);
    }
}
