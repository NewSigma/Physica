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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/LUDecomp.h"

namespace Physica {
    template<Scalar T, int Option, size_t Order = Dynamic>
    class LUSolver {
        using This = LUSolver<T, Option, Order>;
        using MatrixType = DenseMatrix<T, Option, Order, Order>;
        using MatrixB = DenseMatrix<T, MatrixOption::Col | MatrixOption::Vector>;
    private:
        LUDecomp<T, true> lu;
    public:
        LUSolver() = default;
        LUSolver(MatrixType A);
        LUSolver(const This&) = default;
        LUSolver(This&&) noexcept = default;
        ~LUSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void decomp(MatrixType A);
        DenseVector<T, Order> solve(const Vector auto& b) const;
        MatrixB solve(const Matrix auto& b);

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return lu.getOrder(); }
    };

    template<Scalar T, int Option, size_t Order>
    LUSolver<T, Option, Order>::LUSolver(MatrixType A) : lu(std::move(A)) {}

    template<Scalar T, int Option, size_t Order>
    void LUSolver<T, Option, Order>::decomp(MatrixType A) {
        assert(A.isSquare());
        if (A.getRow() != getOrder())
            resize(A.getRow());
        lu.compute(std::move(A));
    }

    template<Scalar T, int Option, size_t Order>
    DenseVector<T, Order> LUSolver<T, Option, Order>::solve(const Vector auto& b) const {
        const size_t order = getOrder();
        assert(b.getLength() == order);
        DenseVector<T, Order> result = lu.getPerm() * b;
        MatrixType working = lu.getMatrixLU();
        for (size_t i = 0; i < order - 1; ++i) {
            auto bottom = working.bottomRows(i + 1);
            result.tail(i + 1) -= bottom.col(i) * result[i];
        }

        for (size_t i = order - 1; i > 0; --i) {
            result[i] /= working(i, i);
            auto top = working.topRows(i);
            result.head(i) -= top.col(i) * result[i];
        }
        result[0] /= working(0, 0);
        return result;
    }

    template<Scalar T, int Option, size_t Order>
    auto LUSolver<T, Option, Order>::solve(const Matrix auto& b) -> MatrixB {
        assert(b.getRow() == getOrder());
        MatrixB result(getOrder(), b.getCol());
        for (size_t i = 0; i < b.getCol(); ++i)
            result.asArray() = solve(result.col(i));
        return result;
    }

    template<Scalar T, int Option, size_t Order>
    void LUSolver<T, Option, Order>::resize(size_t size) {
        lu.resize(size);
    }

    template<Scalar T, int Option, size_t Order>
    void LUSolver<T, Option, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lu.swap(obj.lu);
    }
}
