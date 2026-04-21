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

#include "../RValueMatrix.h"

namespace Physica {
    template<class M, bool ReduceCol>
    class MatrixSum : public RValueVector<MatrixSum<M, ReduceCol>> {
        using This = MatrixSum<M, ReduceCol>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    private:
        const M& mat;
    public:
        MatrixSum(const M& mat_) : mat(mat_) {}
        MatrixSum(const This&) = default;
        MatrixSum(This&&) noexcept = default;
        ~MatrixSum() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        void reverse(const Vector auto& grad) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
    };

    template<class M, bool ReduceCol>
    auto MatrixSum<M, ReduceCol>::calc(size_t index) const -> T {
        if constexpr (ReduceCol)
            return mat.row(index).sum();
        else
            return mat.col(index).sum();
    }

    template<class M, bool ReduceCol>
    auto MatrixSum<M, ReduceCol>::calc_value(size_t index) const -> Tv {
        if constexpr (ReduceCol)
            return mat.values().row(index).sum();
        else
            return mat.values().col(index).sum();
    }

    template<class M, bool ReduceCol>
    void MatrixSum<M, ReduceCol>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        if constexpr (ReduceCol)
            mat.reverse(grad);
        else
            mat.transpose().reverse(grad);
    }

    template<class M, bool ReduceCol>
    size_t MatrixSum<M, ReduceCol>::getLength() const noexcept {
        if constexpr (ReduceCol)
            return mat.getRow();
        else
            return mat.getCol();
    }
}

namespace Physica {
    template<Matrix M, bool ReduceCol>
    class Traits<MatrixSum<M, ReduceCol>> {
    public:
        using ScalarType = M::ScalarType;
        constexpr static size_t SizeAtCompile = ReduceCol ? M::RowAtCompile : M::ColAtCompile;
    };
}
