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
    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    class DenseMatrix;

    template<Matrix M1, Matrix M2>
    class GEMM : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
        constexpr static int ThresholdMKL = 32; // Based on benchmark

        enum InnerDim : char {
            M, K, N
        };
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    private:
        decay_rvalue_t<M1> mat1;
        decay_rvalue_t<M2> mat2;
    public:
        GEMM(M1&& mat1_, M2&& mat2_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        using Base::operator*;
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const;
        void assign_mkl(Matrix auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Matrix auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add(Matrix auto&& target) const noexcept;
        [[nodiscard]] auto compute() const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;

        [[nodiscard]] T trace() const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    private:
        void assign_add_blocking(Matrix auto&& target) const noexcept;
        void assign_add_m(Matrix auto&& target) const noexcept;
        void assign_add_k(Matrix auto&& target) const noexcept;
        void assign_add_n(Matrix auto&& target) const noexcept;
        /* Static members */
        [[nodiscard]] consteval static InnerDim getInnerDim(Matrix auto&& target) noexcept;
        [[nodiscard]] constexpr static size_t bisection(size_t size) noexcept;
        [[nodiscard]] consteval static bool isStaticGEMM() noexcept;
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
        static_assert(getRowAtCompile() != 1 && getColAtCompile() != 1, "[Error]: We reject it and expect higher layers to handle it");
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::operator*(this auto&& self, Vector auto&& v) noexcept {
        using Self = decltype(self);
        assert(self.getCol() == v.getLength());
        return std::forward<Self>(self).getLHS() * (std::forward<Self>(self).getRHS() * std::forward<decltype(v)>(v));
    }

    template<Matrix M1, Matrix M2>
    template<ExecutePolicy P>
    void GEMM<M1, M2>::assign(Matrix auto&& target) const {
        target.assert_assign(*this);
        if constexpr (UseMKL(target)) {
            if (getLHS().getSize() > ThresholdMKL && getRHS().getSize() > ThresholdMKL)
                assign_mkl(target);
            else
                assign_base(target);
        }
        else
            assign_base(target);
    }

    template<Matrix M1, Matrix M2>
    template<ExecutePolicy P>
    void GEMM<M1, M2>::assign_base(Matrix auto&& target) const noexcept {
        if constexpr (isStaticGEMM())
            Base::template assign_base<P>(target);
        else {
            target.zeros();
            assign_add(target);
        }
    }

    template<Matrix M1, Matrix M2>
    template<ExecutePolicy P>
    void GEMM<M1, M2>::assign_add(Matrix auto&& target) const noexcept {
        if constexpr (isStaticGEMM())
            Base::assign_add(target);
        else
            assign_add_blocking(target);
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::compute() const {
        constexpr int Major1 = MatrixMajor::isRowMatrix<M1>() ? MatrixMajor::Col : MatrixMajor::Row;
        constexpr int Major2 = MatrixMajor::isSameMajor<M1, M2>() ? Major1 : MatrixMajor::BothMajor;
        constexpr int Major = Major2 == MatrixMajor::BothMajor ? MatrixMajor::Col : Major2;
        using RtnTy = DenseMatrix<T, Major, Base::getRowAtCompile(), Base::getColAtCompile(), HostAllocator<T>>;
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
    void GEMM<M1, M2>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        if constexpr (ReverseDiff<M1>)
            mat1.reverse(grad * mat2.transpose());
        if constexpr (ReverseDiff<M2>)
            mat2.reverse(mat1.transpose() * grad);
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::trace() const noexcept -> T {
        assert(Base::isSquare());
        if constexpr (canonicalized<M1, M2>()) // Heuristic
            return hadamard(getLHS(), getRHS().transpose()).sum();
        else
            return (getRHS() * getLHS()).trace();
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M1, Matrix M2>
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.mat1);
    }

    template<Matrix M1, Matrix M2>
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.mat2);
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ consteval size_t GEMM<M1, M2>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M1>::getRowAtCompile();
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ consteval size_t GEMM<M1, M2>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M2>::getColAtCompile();
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ consteval int GEMM<M1, M2>::getMajor() noexcept {
        return MatrixMajor::isSameMajor<M1, M2>() ? MatrixMajor::getMajor<M1>() : MatrixMajor::BothMajor;
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_add_blocking(Matrix auto&& target) const noexcept {
        auto blockingM = [&](int size) noexcept -> bool {
            size_t m = getRow();
            bool success = m > size;
            if (success) {
                size_t half = bisection(m);
                (getLHS().topRows(half) * getRHS()).assign_add(target.topRows(half));
                (getLHS().bottomRows(half) * getRHS()).assign_add(target.bottomRows(half));
            }
            return success;
        };

        auto blockingK = [&](int size) noexcept -> bool {
            size_t k = getLHS().getCol();
            bool success = k > size;
            if (success) {
                size_t half = bisection(k);
                (getLHS().leftCols(half) * getRHS().topRows(half)).assign_add(target);
                (getLHS().rightCols(half) * getRHS().bottomRows(half)).assign_add(target);
            }
            return success;
        };

        auto blockingN = [&](int size) noexcept -> bool {
            size_t n = getCol();
            bool success = n > size;
            if (success) {
                size_t half = bisection(n);
                (getLHS() * getRHS().leftCols(half)).assign_add(target.leftCols(half));
                (getLHS() * getRHS().rightCols(half)).assign_add(target.rightCols(half));
            }
            return success;
        };

        constexpr static int BlockingL1 = Base::calcBlockingSize(HostDevAttr::CacheSizeL1D);
        switch (getInnerDim(target)) {
        case InnerDim::M:
            if (blockingK(BlockingL1))
                return;
            return assign_add_m(target);
        case InnerDim::K:
            if constexpr (MatrixMajor::isRowMatrix<M1>() && MatrixMajor::isColMatrix<M2>()) {
                if (blockingM(BlockingL1) || blockingN(BlockingL1))
                    return;
            }
            else {
                if (blockingK(BlockingL1))
                    return;
            }
            return assign_add_k(target);
        case InnerDim::N:
            if (blockingK(BlockingL1))
                return;
            return assign_add_n(target);
        default:
            unreachable();
        }
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_add_m(Matrix auto&& target) const noexcept {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        const size_t k = lhs.getCol();
        const size_t n = getCol();
        for (size_t c = 0; c < n; ++c) {
            auto col = target.col(c);
            for (size_t i = 0; i < k; ++i)
                col += lhs.col(i) * rhs.calc(i, c);
        }
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_add_k(Matrix auto&& target) const noexcept {
        Base::assign_add(target);
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_add_n(Matrix auto&& target) const noexcept {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        const size_t m = getRow();
        const size_t k = lhs.getCol();
        for (size_t r = 0; r < m; ++r) {
            auto row = target.row(r);
            for (size_t i = 0; i < k; ++i)
                row += lhs.calc(r, i) * rhs.row(i);
        }
    }

    template<Matrix M1, Matrix M2>
    consteval auto GEMM<M1, M2>::getInnerDim(Matrix auto&& target) noexcept -> InnerDim {
        using M = decltype(target);
        if constexpr (MatrixMajor::isSameMajor<M, M1>()) {
            if constexpr (MatrixMajor::isSameMajor<M, M2>()) {
                if constexpr (MatrixMajor::isColMatrix<M>())
                    return InnerDim::M;
                else
                    return InnerDim::N;
            }
            else {
                if constexpr (MatrixMajor::isColMatrix<M>())
                    return InnerDim::M;
                else
                    return InnerDim::K;
            }
        }
        else {
            if constexpr (MatrixMajor::isSameMajor<M, M2>()) {
                if constexpr (MatrixMajor::isColMatrix<M>())
                    return InnerDim::K;
                else
                    return InnerDim::N;
            }
            else
                return InnerDim::K;
        }
    }

    template<Matrix M1, Matrix M2>
    constexpr size_t GEMM<M1, M2>::bisection(size_t size) noexcept {
        constexpr int PacketSize = BestPacket<T, Dynamic>::Size;
        size_t result = (size / 2 + (PacketSize - 1)) / PacketSize * PacketSize;
        assert(result < size && "[Error]: Input size is too small");
        return result;
    }

    template<Matrix M1, Matrix M2>
    consteval bool GEMM<M1, M2>::isStaticGEMM() noexcept {
        return std::remove_cvref_t<M1>::getSizeAtCompile() != Dynamic
            && std::remove_cvref_t<M2>::getSizeAtCompile() != Dynamic;
    }

    template<Matrix M1, Matrix M2>
    consteval bool GEMM<M1, M2>::UseMKL(const Matrix auto& target) noexcept {
        using M = decltype(target);
        using T1 = std::remove_cvref_t<M1>;
        using T2 = std::remove_cvref_t<M2>;
        constexpr bool Large1 = T1::getSizeAtCompile() == Dynamic || T1::getSizeAtCompile() > ThresholdMKL;
        constexpr bool Large2 = T2::getSizeAtCompile() == Dynamic || T2::getSizeAtCompile() > ThresholdMKL;
        constexpr bool UseMKL1 = Large1 && Large2;
        constexpr bool UseMKL2 = MatOrTransUseMKL<M, T1>();
        constexpr bool UseMKL3 = MatOrTransUseMKL<M, T2>();
        return UseMKL1 && UseMKL2 && UseMKL3;
    }

    template<Matrix M1, Matrix M2>
    template<Matrix Target, Matrix MaybeTrans>
    consteval bool GEMM<M1, M2>::MatOrTransUseMKL() noexcept {
        if constexpr (instanceof<MaybeTrans, Transpose>)
            return HasMKL() && Internal::EnableLAPACK<Target, decltype(std::declval<MaybeTrans>().transpose())>::value;
        else
            return HasMKL() && Internal::EnableLAPACK<Target, MaybeTrans>::value;
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<GEMM<M1, M2>> {
        using T1 = std::remove_cvref_t<M1>;
        using T2 = std::remove_cvref_t<M2>;
        static_assert(T1::getColAtCompile() == T2::getRowAtCompile() || T1::getColAtCompile() == Dynamic || T2::getRowAtCompile() == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
    };
}

#ifdef PHYSICA_MKL
    #include "GEMM_MKL.h"
#endif
