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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/PLUDecomposition.h"

namespace Physica::Core {
    template<Scalar T, int Option, size_t Order = Dynamic>
    class LUSolver {
        using LUType = PLUDecomposition<T, Option, Order, Order>;
        using MatrixType = LUType::MatrixType;
    private:
        LUType lu;
    public:
        LUSolver() = default;
        LUSolver(MatrixType A);
        LUSolver(const LUSolver&) = default;
        LUSolver(LUSolver&&) noexcept = default;
        ~LUSolver() = default;
        /* Operators */
        LUSolver& operator=(LUSolver obj) noexcept;
        /* Operations */
        void decomposition(MatrixType A);
        template<Vector V>
        DenseVector<T, Order> solve(const V& b) const;
        void swap(LUSolver& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const { return lu.getMatrix().getRow(); }
    };

    template<Scalar T, int Option, size_t Order>
    LUSolver<T, Option, Order>::LUSolver(MatrixType A) : lu(std::move(A)) {}

    template<Scalar T, int Option, size_t Order>
    LUSolver<T, Option, Order>&
    LUSolver<T, Option, Order>::operator=(LUSolver<T, Option, Order> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T, int Option, size_t Order>
    void LUSolver<T, Option, Order>::decomposition(MatrixType A) {
        lu.compute(std::move(A));
    }

    template<Scalar T, int Option, size_t Order>
    template<Vector V>
    DenseVector<T, Order> LUSolver<T, Option, Order>::solve(const V& b) const {
        const size_t order = getOrder();
        DenseVector<T, Order> result(order);
        for (size_t i = 0; i < order; ++i)
            result[i] = b.calc(lu.getBiasOrder()[i]);

        MatrixType working = lu.getMatrix();
        for (size_t i = 0; i < order - 1; ++i) {
            auto bottom = working.bottomRows(i + 1);
            auto tail = result.tail(i + 1);
            tail -= bottom.col(i) * result[i];
        }
        for (size_t i = order - 1; i > 0; --i) {
            result[i] /= working(i, i);
            auto top = working.topRows(i);
            auto head = result.head(i);
            head -= top.col(i) * result[i];
        }
        result[0] /= working(0, 0);
        return result;
    }

    template<Scalar T, int Option, size_t Order>
    void LUSolver<T, Option, Order>::swap(LUSolver& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lu.swap(obj.lu);
    }
}
