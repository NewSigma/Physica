/*
 * Copyright 2023-2025 Weibo He.
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

namespace Physica {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006:56-62
     */
    template<Scalar T, size_t Dim>
    class LineSearch {
    public:
        using VectorType = DenseVector<T, Dim>;

        T maxStepSize;
        T decreaseCondNum;
        T curvatureCondNum;
        size_t maxIteration;
    public:
        LineSearch(T maxStepSize_, T decreaseCondNum_ = 1E-4, T curvatureCondNum_ = 0.1, size_t maxIteration_ = 100);
        LineSearch(const LineSearch&) = default;
        LineSearch(LineSearch&&) noexcept = default;
        ~LineSearch() = default;
        /* Operators */
        LineSearch& operator=(LineSearch obj) noexcept;
        /* Operations */
        [[nodiscard]] T run(
                std::invocable<VectorType> auto fn,
                std::invocable<VectorType> auto grad,
                const VectorType& x,
                const VectorType& gradient,
                const VectorType& direction) const;
        void swap(LineSearch& __restrict obj) noexcept;
    private:
        [[nodiscard]] T zoom(
                std::invocable<VectorType> auto fn,
                std::invocable<VectorType> auto grad,
                const VectorType& x,
                const VectorType& gradient,
                const VectorType& direction,
                T step1,
                T step2) const;
        [[nodiscard]] T interpolate(
                std::invocable<VectorType> auto fn,
                std::invocable<VectorType> auto grad,
                const VectorType& x,
                const VectorType& direction,
                const T& step1,
                const T& step2) const;
    };

    template<Scalar T, size_t Dim>
    LineSearch<T, Dim>::LineSearch(
                T maxStepSize_,
                T decreaseCondNum_,
                T curvatureCondNum_,
                size_t maxIteration_)
            : maxStepSize(maxStepSize_)
            , decreaseCondNum(decreaseCondNum_)
            , curvatureCondNum(curvatureCondNum_)
            , maxIteration(maxIteration_) {
        assert(maxStepSize.isPositive());
        assert(decreaseCondNum.isPositive());
        assert(decreaseCondNum < curvatureCondNum);
        assert(curvatureCondNum < T(1));
    }

    template<Scalar T, size_t Dim>
    LineSearch<T, Dim>& LineSearch<T, Dim>::operator=(LineSearch<T, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T, size_t Dim>
    T LineSearch<T, Dim>::run(
            std::invocable<VectorType> auto fn,
            std::invocable<VectorType> auto grad,
            const VectorType& x,
            const VectorType& gradient,
            const VectorType& direction) const {
        if (gradient.squaredNorm().isSubNormal())
            return 0;
        T step_lower = T(0);
        const T step_upper = maxStepSize;
        T step = maxStepSize / T(2);

        const T phi_0 = fn(x);
        const T diff_phi_0 = gradient * direction;
        assert(diff_phi_0.isNegative());

        T last_y = phi_0;
        VectorType x1 = VectorType(x.getLength());
        size_t i = 0;
        while (true) {
            x1 = x + step * direction;
            const T y = fn(x1);
            const bool violatesWolfe = y > (phi_0 + decreaseCondNum * step * diff_phi_0);
            const bool isIncreased = (y >= last_y) && (i > 0);
            if (violatesWolfe || isIncreased)
                return zoom(fn, grad, x, gradient, direction, step_lower, step);

            const T diff_phi = grad(x1) * direction;
            if (abs(diff_phi) <= -curvatureCondNum * diff_phi_0)
                return step;
            if (!diff_phi.isNegative())
                return zoom(fn, grad, x, gradient, direction, step, step_lower);
            
            step_lower = step;
            step = (step_lower + step_upper) / T(2);
            if (abs(step_upper - step) < std::numeric_limits<T>::epsilon())
                return step_upper;
            if (++i > maxIteration)
                throw BadConvergenceException("Exceed max iteration of LineSearch");
        }
    }

    template<Scalar T, size_t Dim>
    void LineSearch<T, Dim>::swap(LineSearch& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        maxStepSize.swap(obj.maxStepSize);
        decreaseCondNum.swap(obj.decreaseCondNum);
        curvatureCondNum.swap(obj.curvatureCondNum);
        std::swap(maxIteration, obj.maxIteration);
    }
    /**
     * \param step1 does not necessarily less than \param step2
     */
    template<Scalar T, size_t Dim>
    T LineSearch<T, Dim>::zoom(
            std::invocable<VectorType> auto fn,
            std::invocable<VectorType> auto grad,
            const VectorType& x,
            const VectorType& gradient,
            const VectorType& direction,
            T step1,
            T step2) const {
        const T phi_0 = fn(x);
        const T diff_phi_0 = gradient * direction;

        T last_y = fn(VectorType(x + step1 * direction));
        VectorType x1 = VectorType(x.getLength());
        while (true) {
            const T step = interpolate(fn, grad, x, direction, step1, step2);

            x1 = x + step * direction;
            const T y = fn(x1);
            const bool violatesWolfe = y > (phi_0 + decreaseCondNum * step * diff_phi_0);
            const bool isIncreased = y >= last_y;
            if (violatesWolfe || isIncreased)
                step2 = step;
            else {
                const T diff_phi = grad(x1) * direction;
                if (abs(diff_phi) <= -curvatureCondNum * diff_phi_0)
                    return step;
                if (!(diff_phi * (step2 - step1)).isNegative())
                    step2 = step1;
                step1 = step;
                last_y = y;
            }

            if (abs(step2 - step1) < std::numeric_limits<T>::epsilon())
                return step1;
        }
    }

    template<Scalar T, size_t Dim>
    T LineSearch<T, Dim>::interpolate(
            std::invocable<VectorType> auto fn,
            std::invocable<VectorType> auto grad,
            const VectorType& x,
            const VectorType& direction,
            const T& step1,
            const T& step2) const {
        assert(step1 != step2);
        const VectorType x1 = x + step1 * direction;
        const VectorType x2 = x + step2 * direction;
        const T diff1 = grad(x1) * direction;
        const T diff2 = grad(x2) * direction;
        const T delta_step = step1 - step2;
        const T d1 = diff1 + diff2 - T(3) * (fn(x1) - fn(x2)) / delta_step;
        const T squared_d2 = square(d1) - diff1 * diff2;
        assert(!squared_d2.isNegative());
        T d2 = sqrt(squared_d2);
        if (delta_step.isNegative())
            d2 = -d2;
        return step1 - delta_step * (diff1 + d2 - d1) / (diff1 - diff2 + T(2) * d2);
    }
}
