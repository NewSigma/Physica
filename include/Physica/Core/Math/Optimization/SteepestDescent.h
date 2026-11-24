/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        void init(VectorType initial, std::invocable<VectorType> auto fn, std::invocable<VectorType> auto grad);
        void step(std::invocable<VectorType> auto fn, std::invocable<VectorType> auto grad);
        void swap(SteepestDescent& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept;
        [[nodiscard]] const VectorType& getGradG() const noexcept { return gradG; }
        [[nodiscard]] const VectorType& getArgX() const noexcept { return nowX; }
        [[nodiscard]] T getObjectiveValue() const noexcept { return nowY; }
    };
}

#include "OptimizationImpl/SteepestDescentImpl.h"
