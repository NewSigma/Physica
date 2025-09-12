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

#include "../RValueMatrixImpl.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Matrix M>
    class TrigUpper<M> : public RValueMatrix<TrigUpper<M>> {
        using This = TrigUpper<M>;
        using Base = RValueMatrix<This>;

        M& mat;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Trv;
        using typename Base::Tv;
    public:
        TrigUpper(M& mat_);
        TrigUpper(const This&) = default;
        TrigUpper(This&&) noexcept = default;
        ~TrigUpper() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] T det() const;
        [[nodiscard]] Tr lnAbsDet() const;
        [[nodiscard]] Trv sgndet() const;

        [[nodiscard]] VectorND<T> toDT();
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getRow() const noexcept { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getCol(); }
    };

    template<Matrix M>
    TrigUpper<M>::TrigUpper(M& mat_) : mat(mat_) {}

    template<Matrix M>
    auto TrigUpper<M>::calc(size_t row, size_t col) const -> T {
        if (row > col)
            return T(0);
        return mat.calc(row, col);
    }

    template<Matrix M>
    auto TrigUpper<M>::calc_value(size_t row, size_t col) const -> Tv {
        if (row > col)
            return Tv(0);
        return mat.calc_value(row, col);
    }

    template<Matrix M>
    auto TrigUpper<M>::det() const -> T {
        return Base::diag().prod();
    }

    template<Matrix M>
    auto TrigUpper<M>::lnAbsDet() const -> Tr {
        assert(!Base::diag().prod().isZero() && "[Error]: Singular matrix");
        return ln(abs(Base::diag())).sum();
    }

    template<Matrix M>
    auto TrigUpper<M>::sgndet() const -> Trv {
        static_assert(!T::isComplex, "[Error]: sgndet() is not well defined for complex matrix");
        return unit(mat.diag()).prod();
    }

    template<Matrix M>
    auto TrigUpper<M>::toDT() -> VectorND<T> {
        const size_t length = getRow();
        VectorND<T> diag(length);
        for (size_t i = 0; i < length; ++i) {
            diag[i] = calc(i, i);
            mat.row(i).tail(i) *= reciprocal(diag[i]);
        }
        return diag;
    }

    template<Matrix M>
    class TrigLower<M> : public RValueMatrix<TrigLower<M>> {
        using This = TrigLower<M>;
        using Base = RValueMatrix<This>;

        M& mat;
    public:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Trv;
        using typename Base::Tv;
    public:
        TrigLower(M& mat_);
        TrigLower(const This&) = default;
        TrigLower(This&&) noexcept = default;
        ~TrigLower() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;
        [[nodiscard]] Trv sgndet() const;

        [[nodiscard]] T det() const;
        [[nodiscard]] Tr lnAbsDet() const;
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getRow() const noexcept { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getCol(); }
    };

    template<Matrix M>
    TrigLower<M>::TrigLower(M& mat_) : mat(mat_) {
        assert(mat.isSquare());
    }

    template<Matrix M>
    auto TrigLower<M>::calc(size_t row, size_t col) const -> T {
        if (row < col)
            return T(0);
        return mat.calc(row, col);
    }

    template<Matrix M>
    auto TrigLower<M>::calc_value(size_t row, size_t col) const -> Tv {
        if (row < col)
            return Tv(0);
        return mat.calc_value(row, col);
    }

    template<Matrix M>
    auto TrigLower<M>::sgndet() const -> Trv {
        return unit(mat.diag()).prod();
    }

    template<Matrix M>
    auto TrigLower<M>::det() const -> T {
        return Base::diag().prod();
    }

    template<Matrix M>
    auto TrigLower<M>::lnAbsDet() const -> Tr {
        assert(!Base::diag().prod().isZero() && "[Error]: Singular matrix");
        return ln(abs(Base::diag())).sum();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<TrigUpper<M>> : public Traits<M> {
    public:
        constexpr static int Option = MatrixOption::getMajor<M>() | MatrixOption::AnyStorage;
    };

    template<Matrix M>
    class Traits<TrigLower<M>> : public Traits<M> {
    public:
        constexpr static int Option = MatrixOption::getMajor<M>() | MatrixOption::AnyStorage;
    };
}

#include "TrigGEMM.h"
#include "InvTrigGEMM.h"
#include "TrigDiagMM.h"
