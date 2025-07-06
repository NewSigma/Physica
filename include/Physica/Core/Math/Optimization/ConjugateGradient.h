/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "OptimizationImpl/LineSearch.h"

namespace Physica {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006:121-122
     */
    template<Scalar T, size_t Dim>
    class ConjugateGradient {
        using This = ConjugateGradient<T, Dim>;
        using VectorType = DenseVector<T, Dim>;

        VectorType gradG;
        VectorType direction;
        VectorType nowX;
        T squaredGradNorm;
        LineSearch<T, Dim> lineSearch;
        size_t iteration;
    public:
        ConjugateGradient(T maxStepSize);
        ConjugateGradient(const This&) = default;
        ConjugateGradient(This&&) noexcept = default;
        ~ConjugateGradient() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void init(VectorType initial, std::invocable<VectorType> auto func, std::invocable<VectorType> auto grad);
        void step(std::invocable<VectorType> auto func, std::invocable<VectorType> auto grad);
         T solve(T epsilon, VectorType initial, std::invocable<VectorType> auto func, std::invocable<VectorType> auto grad);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] inline size_t getDim() const noexcept;
        [[nodiscard]] const VectorType& getGradG() const noexcept { return gradG; }
        [[nodiscard]] const VectorType& getArgX() const noexcept { return nowX; }
    };

    template<Scalar T, size_t Dim>
    ConjugateGradient<T, Dim>::ConjugateGradient(T maxStepSize)
            : lineSearch(maxStepSize)
            , iteration(0) {}

    template<Scalar T, size_t Dim>
    void ConjugateGradient<T, Dim>::init(VectorType initial, std::invocable<VectorType> auto, std::invocable<VectorType> auto grad) {
        nowX = std::move(initial);
        gradG = grad(nowX);
        direction = -gradG;
        squaredGradNorm = gradG.squaredNorm();
        iteration = 0;
    }

    template<Scalar T, size_t Dim>
    void ConjugateGradient<T, Dim>::step(std::invocable<VectorType> auto func, std::invocable<VectorType> auto grad) {
        const size_t dim = nowX.getLength();
        const T stepSize = lineSearch.run(func, grad, nowX, gradG, direction);
        nowX += stepSize * direction;

        const VectorType gradG1 = grad(nowX);
        const T squaredGradNorm1 = gradG1.squaredNorm();
        if (++iteration == dim) {
            direction = -gradG1;
            iteration = 0;
        }
        else {
            const T alpha = (squaredGradNorm1 - gradG * gradG1) / squaredGradNorm; //Polak-Ribière method[1]
            direction = -gradG1 + alpha * direction;
        }
        gradG = gradG1;
        squaredGradNorm = squaredGradNorm1;
    }

    template<Scalar T, size_t Dim>
    T ConjugateGradient<T, Dim>::solve(T epsilon, VectorType initial, std::invocable<VectorType> auto func, std::invocable<VectorType> auto grad) {
        init(std::move(initial), func, grad);
        while (squaredGradNorm > epsilon) {
            step(func, grad);
        }
        return func(nowX);
    }

    template<Scalar T, size_t Dim>
    void ConjugateGradient<T, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        gradG.swap(obj.gradG);
        direction.swap(obj.direction);
        nowX.swap(obj.nowX);
        squaredGradNorm.swap(obj.squaredGradNorm);
        lineSearch.swap(obj.lineSearch);
        std::swap(iteration, obj.iteration);
    }

    template<Scalar T, size_t Dim>
    inline size_t ConjugateGradient<T, Dim>::getDim() const noexcept {
        if constexpr (Dim != Dynamic)
            return Dim;
        return gradG.getLength();
    }
}
