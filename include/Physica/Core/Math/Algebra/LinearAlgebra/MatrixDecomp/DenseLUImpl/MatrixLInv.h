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

#include "MatrixL.h"

namespace Physica {
    template<Scalar T>
    class Inverse<LUMatrixL<T>> : public RValueMatrix<Inverse<LUMatrixL<T>>> {
        using This = Inverse<LUMatrixL<T>>;
        using Base = RValueMatrix<This>;
        using typename Base::Trv;

        const LUMatrixL<T>& matL;
    public:
        Inverse(const LUMatrixL<T>& matL_);
        Inverse(const This&) = delete;
        Inverse(This&&) noexcept = delete;
        ~Inverse() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return matL.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return matL.getCol(); }
    };

    template<Scalar T>
    Inverse<LUMatrixL<T>>::Inverse(const LUMatrixL<T>& matL_) : matL(matL_) {}

    template<Scalar T>
    void Inverse<LUMatrixL<T>>::assign(Matrix auto& target) const {
        target.assert_assign(*this);
        target = -matL;
        target.diag() = Trv(1);
        for (size_t i = 1; i < getCol() - 1; ++i) {
            auto corner = target.bottomLeftCorner(i + 1, i);

            auto row = target.row(i);
            auto col = target.col(i);
            auto head = row.head(i);
            auto tail = col.tail(i + 1);
            corner += tail * head.transpose();
        }
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<Inverse<LUMatrixL<T>>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::Col;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
