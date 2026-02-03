/*
 * Copyright 2024-2026 Weibo He.
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
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    class DenseMatrix;

    template<Matrix M1, Matrix M2>
    class GEMM : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    public:
        constexpr static int Critical = 32;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    private:
        LazyDestroy<M1> mat1;
        LazyDestroy<M2> mat2;
    public:
        GEMM(M1&& mat1_, M2&& mat2_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const noexcept;
        [[nodiscard]] auto compute() const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    private:
        /* Static members */
        [[nodiscard]] consteval static bool UseMKL(const Matrix auto& target) noexcept;
        template<Matrix Target, Matrix MaybeTrans>
        consteval static bool MatOrTransUseMKL() noexcept;
        /* Friends */
        friend class device_obj<This>;
    };

    template<Matrix M1, Matrix M2>
    GEMM<M1, M2>::GEMM(M1&& mat1_, M2&& mat2_) : mat1(std::forward<M1>(mat1_)), mat2(std::forward<M2>(mat2_)) {
        assert(mat1.getRow() > 0);
        assert(mat1.getCol() > 0);
        assert(mat2.getCol() > 0);
        assert(mat1.getCol() == mat2.getRow());
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        target.assert_assign(*this);
        if constexpr (UseMKL(target)) {
            if (getLHS().getSize() > Critical && getRHS().getSize() > Critical)
                assign_mkl(target);
            else
                Base::assign_base(target);
        }
        else
            Base::assign_base(target);
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::compute() const {
        constexpr static bool SameMajor = MatrixMajor::isSameMajor<M1, M2>();
        constexpr static bool RowMajor = MatrixMajor::isRowMatrix<M1>();
        constexpr static auto Major1 = RowMajor ? MatrixMajor::Col : MatrixMajor::Row;
        constexpr static auto Major = SameMajor ? Major1 : MatrixMajor::BothMajor;
        constexpr static int Option = Major == MatrixMajor::BothMajor ? MatrixMajor::Col : Major;
        using RtnTy = DenseMatrix<T, Option, Base::RowAtCompile, Base::ColAtCompile, HostAllocator<T>>;
        return RtnTy(*this);
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        if constexpr (IsIntelLLVM()) {
            T result(0);
            for (size_t i = 0; i < mat1.getCol(); ++i)
                result += mat1.calc(row, i) * mat2.calc(i, col);
            return result;
        }
        else // Intel LLVM does not optimize the following code well
            return mat1.row(row) * mat2.col(col);
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        return (mat1.values() * mat2.values()).calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        if constexpr (ReverseDiff<M1>)
            mat1.reverse(grad * mat2.transpose());
        if constexpr (ReverseDiff<M2>)
            mat2.reverse(mat1.transpose() * grad);
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::values() const noexcept {
        return mat1.values() * mat2.values();
    }

    template<Matrix M1, Matrix M2>
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.mat1);
    }

    template<Matrix M1, Matrix M2>
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.mat2);
    }

    template<Matrix M1, Matrix M2>
    consteval bool GEMM<M1, M2>::UseMKL(const Matrix auto& target) noexcept {
        using M = decltype(target);
        using T1 = std::remove_cvref_t<M1>;
        using T2 = std::remove_cvref_t<M2>;
        constexpr bool Large1 = T1::SizeAtCompile == Dynamic || T1::SizeAtCompile > Critical;
        constexpr bool Large2 = T2::SizeAtCompile == Dynamic || T2::SizeAtCompile > Critical;
        constexpr bool UseMKL1 = Large1 && Large2;
        constexpr bool UseMKL2 = MatOrTransUseMKL<M, T1>();
        constexpr bool UseMKL3 = MatOrTransUseMKL<M, T2>();
        return UseMKL1 && UseMKL2 && UseMKL3;
    }

    template<Matrix M1, Matrix M2>
    template<Matrix Target, Matrix MaybeTrans>
    consteval bool GEMM<M1, M2>::MatOrTransUseMKL() noexcept {
        if constexpr (instanceof<Transpose, MaybeTrans>)
            return Internal::EnableMKL<Target, decltype(std::declval<MaybeTrans>().getExpr())>::value;
        else
            return Internal::EnableMKL<Target, MaybeTrans>::value;
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<GEMM<M1, M2>> {
        using T1 = std::remove_cvref_t<M1>;
        using T2 = std::remove_cvref_t<M2>;
        static_assert(T1::ColAtCompile == T2::RowAtCompile || T1::ColAtCompile == Dynamic || T2::RowAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        constexpr static int Option = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = T1::RowAtCompile;
        constexpr static size_t ColAtCompile = T2::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#ifdef PHYSICA_MKL
    #include "GEMM_MKL.h"
#endif
