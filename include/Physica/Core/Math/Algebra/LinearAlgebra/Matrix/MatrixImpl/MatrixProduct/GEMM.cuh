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

#include "Physica/PlainStruct.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "GEMM.h"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class device_obj<GEMM<M1, M2>> : public device_obj<RValueMatrix<GEMM<M1, M2>>> {
        using host_obj = GEMM<M1, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref1 = add_device_obj<M1>::type;
        using Ref2 = add_device_obj<M2>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M1>>> mat1;
        PlainStruct<add_device_obj_t<std::remove_reference_t<M2>>> mat2;
    public:
        __host__ __device__ device_obj(Ref1 mat1_, Ref2 mat2_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator*;
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] __host__ __device__ auto operator*(this auto&&, Vector auto&& v) noexcept;
        /* Operations */
        __host__ __device__ void assign(Matrix auto& target) const;
        __host__ __device__ void assign_add(Matrix auto& target) const;

        [[nodiscard]] auto compute() const;
        [[nodiscard]] __device__ T calc(size_t, size_t) const { noImpl("GEMM.calc() is low performance and should be avoided"); }

        void reverse(const Matrix auto& grad) const noexcept;

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getRHS().getCol(); }
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
    private:
        template<bool AssignAdd>
        void assign_impl_cublas(Matrix auto& target) const;
    };

    template<Matrix M1, Matrix M2>
    __host__ __device__ device_obj<GEMM<M1, M2>>::device_obj(Ref1 mat1_, Ref2 mat2_) : mat1(asStruct(mat1_)), mat2(asStruct(mat2_)) {
        assert(mat1_.getCol() == mat2_.getRow());
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto device_obj<GEMM<M1, M2>>::operator*(this auto&& self, Vector auto&& v) noexcept {
        using Self = decltype(self);
        assert(self.getCol() == v.getLength());
        return std::forward<Self>(self).getLHS() * (std::forward<Self>(self).getRHS() * std::forward<decltype(v)>(v));
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ void device_obj<GEMM<M1, M2>>::assign(Matrix auto& target) const {
        if constexpr (IsHost())
            assign_impl_cublas<false>(target);
        else
            noImpl("No device GEMM support");
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ void device_obj<GEMM<M1, M2>>::assign_add(Matrix auto& target) const {
        if constexpr (IsHost())
            assign_impl_cublas<true>(target);
        else
            noImpl("No device GEMM support");
    }

    template<Matrix M1, Matrix M2>
    auto device_obj<GEMM<M1, M2>>::compute() const {
        using RetTy = decltype(std::declval<host_obj>().compute())::device_obj_type;
        return RetTy(*this);
    }

    template<Matrix M1, Matrix M2>
    void device_obj<GEMM<M1, M2>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        if constexpr (ReverseDiff<M1>)
            getLHS().reverse(grad * getRHS().values().transpose());
        if constexpr (ReverseDiff<M2>)
            getRHS().reverse(getLHS().values().transpose() * grad);
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto device_obj<GEMM<M1, M2>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto&& device_obj<GEMM<M1, M2>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref1>(self.mat1.getDerived());
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto&& device_obj<GEMM<M1, M2>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref2>(self.mat2.getDerived());
    }

    template<Matrix M1, Matrix M2>
    template<bool AssignAdd>
    void device_obj<GEMM<M1, M2>>::assign_impl_cublas(Matrix auto& target) const {
        using Tm = decltype(std::declval<T>().toCUDA());
        constexpr bool IsDeviceMatrix = getLHS().getSizeAtCompile() == Dynamic && getRHS().getSizeAtCompile() == Dynamic;
        constexpr bool isTranspose1 = instanceof<Transpose, M1>;
        constexpr bool isTranspose2 = instanceof<Transpose, M2>;
        static_assert(IsDeviceMatrix, "[Error]: Fixed matrix is on host, pass it to device before calling cuBLAS");
        static_assert(MatrixMajor::isColMatrix<M1>() && MatrixMajor::isColMatrix<M2>(), "[Error]: cuBLAS uses column major");
        static_assert(!Diffable<This>);
        target.assert_assign(*this);

        const auto op1 = isTranspose1 ? CUBLAS_OP_T : CUBLAS_OP_N;
        const auto op2 = isTranspose2 ? CUBLAS_OP_T : CUBLAS_OP_N;
        const size_t r = getRow();
        const size_t c = getCol();
        const Tm alpha = 1;
        const Tm* pAlpha = &alpha;

        Tm beta;
        if constexpr (AssignAdd)
            beta = 1;
        else
            beta = 0;
        const Tm* pBeta = &beta;

        const Tm* A{};
        size_t k{};
        if constexpr (isTranspose1) {
            k = getLHS().getRow();
            A = (Tm*)getLHS().getExpr().data();
        }
        else {
            k = getLHS().getCol();
            A = (Tm*)getLHS().data();
        }

        const Tm* B{};
        if constexpr (isTranspose2)
            B = (Tm*)getRHS().getExpr().data();
        else
            B = (Tm*)getRHS().data();

        Tm* C{};
        if constexpr (Diffable<decltype(target)>)
            C = (Tm*)target.data().value_ptr();
        else
            C = (Tm*)target.data();

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (T::Prec == Float16)
            check(cublasHgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        else if constexpr (T::Prec == Float32)
            check(cublasSgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        else {
            static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
            check(cublasDgemm_64(ctx, op1, op2, r, c, k, pAlpha, A, r, B, k, pBeta, C, r));
        }
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<device_obj<GEMM<M1, M2>>> : public Traits<GEMM<M1, M2>> {
        static_assert(!is_device_obj<M1>::value);
        static_assert(!is_device_obj<M2>::value);
    };
}
