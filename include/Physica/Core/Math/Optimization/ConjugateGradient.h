/*
 * Copyright 2021-2023 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "OptimizationImpl/LineSearch.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.121-122
     */
    template<class ScalarType, size_t Dim>
    class ConjugateGradient {
        using VectorType = Vector<ScalarType, Dim>;

        VectorType gradG;
        VectorType direction;
        VectorType nowX;
        ScalarType squaredGradNorm;
        LineSearch<ScalarType, Dim> lineSearch;
        size_t iteration;
    public:
        ConjugateGradient(ScalarType maxStepSize);
        ConjugateGradient(const ConjugateGradient&) = default;
        ConjugateGradient(ConjugateGradient&&) noexcept = default;
        ~ConjugateGradient() = default;
        /* Operators */
        ConjugateGradient& operator=(ConjugateGradient obj) noexcept;
        /* Operations */
        template<class Functor, class GradFunctor> void init(VectorType initial, Functor func, GradFunctor grad);
        template<class Functor, class GradFunctor> void step(Functor func, GradFunctor grad);
        template<class Functor, class GradFunctor> ScalarType solve(ScalarType epsilon, VectorType initial, Functor func, GradFunctor grad);
        void swap(ConjugateGradient& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] inline size_t getDim() const noexcept;
        [[nodiscard]] const VectorType& getGradG() const noexcept { return gradG; }
        [[nodiscard]] const VectorType& getArgX() const noexcept { return nowX; }
    };

    template<class ScalarType, size_t Dim>
    ConjugateGradient<ScalarType, Dim>::ConjugateGradient(ScalarType maxStepSize)
            : lineSearch(maxStepSize)
            , iteration(0) {}

    template<class ScalarType, size_t Dim>
    ConjugateGradient<ScalarType, Dim>&
    ConjugateGradient<ScalarType, Dim>::operator=(ConjugateGradient<ScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    void ConjugateGradient<ScalarType, Dim>::init(VectorType initial, [[maybe_unused]] Functor func, GradFunctor grad) {
        nowX = std::move(initial);
        gradG = grad(nowX);
        direction = -gradG;
        squaredGradNorm = gradG.squaredNorm();
        iteration = 0;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    void ConjugateGradient<ScalarType, Dim>::step(Functor func, GradFunctor grad) {
        const size_t dim = nowX.getLength();
        const ScalarType stepSize = lineSearch.run(func, grad, nowX, gradG, direction);
        nowX += stepSize * direction;

        const VectorType gradG1 = grad(nowX);
        const ScalarType squaredGradNorm1 = gradG1.squaredNorm();
        if (++iteration == dim) {
            direction = -gradG1;
            iteration = 0;
        }
        else {
            const ScalarType alpha = (squaredGradNorm1 - gradG * gradG1) / squaredGradNorm; //Polak-Ribière method[1]
            direction = -gradG1 + alpha * direction;
        }
        gradG = gradG1;
        squaredGradNorm = squaredGradNorm1;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    ScalarType ConjugateGradient<ScalarType, Dim>::solve(ScalarType epsilon, VectorType initial, Functor func, GradFunctor grad) {
        init(std::move(initial), func, grad);
        while (squaredGradNorm > epsilon) {
            step(func, grad);
        }
        return func(nowX);
    }

    template<class ScalarType, size_t Dim>
    void ConjugateGradient<ScalarType, Dim>::swap(ConjugateGradient& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        gradG.swap(obj.gradG);
        direction.swap(obj.direction);
        nowX.swap(obj.nowX);
        squaredGradNorm.swap(obj.squaredGradNorm);
        lineSearch.swap(obj.lineSearch);
        std::swap(iteration, obj.iteration);
    }

    template<class ScalarType, size_t Dim>
    inline size_t ConjugateGradient<ScalarType, Dim>::getDim() const noexcept {
        if constexpr (Dim != Dynamic)
            return Dim;
        return gradG.getLength();
    }
}
