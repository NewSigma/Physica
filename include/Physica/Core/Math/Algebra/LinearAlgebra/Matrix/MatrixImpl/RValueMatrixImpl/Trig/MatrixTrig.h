/*
 * Copyright 2025-2026 Weibo He.
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
    template<Matrix M, bool Upper, bool Unit>
    class MatrixTrig<M, Upper, Unit> : public RValueMatrix<MatrixTrig<M, Upper, Unit>> {
        using This = MatrixTrig<M, Upper, Unit>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
        using typename Base::Trv;
        using typename Base::Tm;
    private:
        LazyDestroy<M&&> mat;
    public:
        MatrixTrig(M mat_);
        MatrixTrig(const This&) = default;
        MatrixTrig(This&&) noexcept = default;
        ~MatrixTrig() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const;
        void assign_mkl(Matrix auto&& target) const;
        // TODO: assign_base() lower it to memcpy

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] auto det() const;
        [[nodiscard]] auto lnAbsDet() const;
        [[nodiscard]] auto sgndet() const;
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getRow() const noexcept { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getCol(); }
    };

    template<Matrix M, bool Upper, bool Unit>
    MatrixTrig<M, Upper, Unit>::MatrixTrig(M mat) : mat(std::forward<M>(mat)) {}

    template<Matrix M, bool Upper, bool Unit>
    template<ExecutePolicy P>
    void MatrixTrig<M, Upper, Unit>::assign(Matrix auto&& target) const {
        Base::template assign_base<P>(target);
    }

    template<Matrix M, bool Upper, bool Unit>
    auto MatrixTrig<M, Upper, Unit>::calc(size_t row, size_t col) const -> T {
        if constexpr (Upper) {
            if (row > col)
                return Trv(0);
        }
        else {
            if (row < col)
                return Trv(0);
        }

        if constexpr (Unit)
            if (row == col)
                return Trv(1);
        return mat.calc(row, col);
    }

    template<Matrix M, bool Upper, bool Unit>
    auto MatrixTrig<M, Upper, Unit>::calc_value(size_t row, size_t col) const -> Tv {
        if constexpr (Unit)
            if (row == col)
                return Trv(1);

        if constexpr (Upper) {
            if (row > col)
                return Trv(0);
        }
        else {
            if (row < col)
                return Trv(0);
        }
        return mat.calc_value(row, col);
    }

    template<Matrix M, bool Upper, bool Unit>
    auto MatrixTrig<M, Upper, Unit>::det() const {
        if constexpr (Unit)
            return Trv(1);
        else
            return Base::diag().prod();
    }

    template<Matrix M, bool Upper, bool Unit>
    auto MatrixTrig<M, Upper, Unit>::lnAbsDet() const {
        if constexpr (Unit)
            return Trv(0);
        else {
            assert(!Base::diag().prod().isZero() && "[Error]: Singular matrix");
            return ln(abs(Base::diag())).sum();
        }
    }

    template<Matrix M, bool Upper, bool Unit>
    auto MatrixTrig<M, Upper, Unit>::sgndet() const {
        if constexpr (Unit)
            return Trv(1);
        else
            return unit(mat.diag()).prod();
    }
}

namespace Physica {
    template<Matrix M, bool Upper_, bool Unit_>
    class Traits<MatrixTrig<M, Upper_, Unit_>> : public Traits<M> {
    public:
        constexpr static int Option = MatrixOption::getMajor<M>();
        constexpr static bool Upper = Upper_;
        constexpr static bool Unit = Unit_;

        using ExprType = std::remove_cvref<M>::type;
    };
}

#ifdef PHYSICA_MKL
    #include "MatrixTrig_MKL.h"
#endif
#include "GEMM.h"
#include "Inverse.h"
#include "InvGEMV.h"
#include "InvGEMM.h"
#include "DiagGEMM.h"
