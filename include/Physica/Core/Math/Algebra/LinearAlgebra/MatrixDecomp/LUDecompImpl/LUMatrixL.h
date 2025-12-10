/*
 * Copyright 2025 Weibo He.
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

#include "../LUDecomp.h"

namespace Physica {
    template<Scalar T>
    class LUMatrixL : public RValueMatrix<LUMatrixL<T>> {
        using This = LUMatrixL<T>;
        using Base = RValueMatrix<This>;
        using typename Base::Trv;

        const DenseMatrix<T>& matrixLU;
    public:
        template<bool Pivot>
        LUMatrixL(const LUDecomp<T, Pivot>& lu);
        LUMatrixL(const This&) = delete;
        LUMatrixL(This&&) noexcept = delete;
        ~LUMatrixL() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return matrixLU.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return matrixLU.getCol(); }
    };

    template<Scalar T>
    template<bool Pivot>
    LUMatrixL<T>::LUMatrixL(const LUDecomp<T, Pivot>& lu) : matrixLU(lu.getMatrixLU()) {}

    template<Scalar T>
    T LUMatrixL<T>::calc(size_t row, size_t col) const {
        if (row == col)
            return T(1);
        return matrixLU.tril().calc(row, col);
    }

    template<Scalar T>
    void LUMatrixL<T>::assign(Matrix auto& target) const {
        target = matrixLU.tril();
        target.diag() = Trv(1);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<LUMatrixL<T>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::Col;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#include "LUMatrixLInv.h"
