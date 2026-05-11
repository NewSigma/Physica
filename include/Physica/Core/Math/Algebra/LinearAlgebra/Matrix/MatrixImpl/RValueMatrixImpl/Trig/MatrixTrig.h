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
    private:
        LazyDestroy<M> mat;
    public:
        MatrixTrig(M&& mat_);
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

        [[nodiscard]] auto det() const;
        [[nodiscard]] auto lnAbsDet() const;
        [[nodiscard]] auto sgndet() const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
    };

    template<Matrix M, bool Upper, bool Unit>
    MatrixTrig<M, Upper, Unit>::MatrixTrig(M&& mat) : mat(std::forward<M>(mat)) {}

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

    template<Matrix M, bool Upper, bool Unit>
    auto MatrixTrig<M, Upper, Unit>::values(this auto&& self) noexcept {
        using Value = decltype(std::forward<decltype(self)>(self).getExpr().values());
        return MatrixTrig<Value, Upper, Unit>(std::forward<decltype(self)>(self).getExpr().values());
    }

    template<Matrix M, bool Upper, bool Unit>
    auto&& MatrixTrig<M, Upper, Unit>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, bool Upper, bool Unit>
    size_t MatrixTrig<M, Upper, Unit>::getRow() const noexcept {
        if constexpr (Upper)
            return std::min(mat.getRow(), mat.getCol());
        else
            return mat.getRow();
    }

    template<Matrix M, bool Upper, bool Unit>
    size_t MatrixTrig<M, Upper, Unit>::getCol() const noexcept {
        if constexpr (Upper)
            return mat.getCol();
        else
            return std::min(mat.getRow(), mat.getCol());
    }

    template<Matrix M, bool Upper, bool Unit>
    __host__ __device__ consteval size_t MatrixTrig<M, Upper, Unit>::getRowAtCompile() noexcept {
        using M1 = std::remove_cvref_t<M>;
        if constexpr (Upper)
            return std::min(M1::getRowAtCompile(), M1::getColAtCompile());
        else
            return M1::getRowAtCompile();
    }

    template<Matrix M, bool Upper, bool Unit>
    __host__ __device__ consteval size_t MatrixTrig<M, Upper, Unit>::getColAtCompile() noexcept {
        using M1 = std::remove_cvref_t<M>;
        if constexpr (Upper)
            return M1::getColAtCompile();
        else
            return std::min(M1::getRowAtCompile(), M1::getColAtCompile());
    }
}

namespace Physica {
    template<Matrix M, bool Upper_, bool Unit_>
    class Traits<MatrixTrig<M, Upper_, Unit_>> : public Traits<M> {
    public:
        constexpr static int Major = MatrixMajor::getMajor<M>();
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
