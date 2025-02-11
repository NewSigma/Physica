/*
 * Copyright 2020-2023 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "OptimizationImpl/LineSearch.h"

namespace Physica {
    template<Scalar T, size_t Dim>
    class SteepestDescent {
    public:
        using VectorType = DenseVector<T, Dim>;
    private:
        VectorType gradG;
        VectorType tryX;
        VectorType nowX;
        T nowY;
        LineSearch<T, Dim> lineSearch;
    public:
        SteepestDescent(T maxStepSize, T decreaseCondNum_, T curvatureCondNum_);
        SteepestDescent(const SteepestDescent&) = default;
        SteepestDescent(SteepestDescent&&) noexcept = default;
        ~SteepestDescent() = default;
        /* Operators */
        SteepestDescent& operator=(SteepestDescent obj) noexcept;
        /* Operations */
        template<class Functor, class GradFunctor> void init(VectorType initial, Functor func, GradFunctor grad);
        template<class Functor, class GradFunctor> void step(Functor func, GradFunctor grad);
        void swap(SteepestDescent& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] inline size_t getDim() const noexcept;
        [[nodiscard]] const VectorType& getGradG() const noexcept { return gradG; }
        [[nodiscard]] const VectorType& getArgX() const noexcept { return nowX; }
        [[nodiscard]] T getObjectiveValue() const noexcept { return nowY; }
    };
}

#include "OptimizationImpl/SteepestDescentImpl.h"
