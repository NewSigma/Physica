/*
 * Copyright 2026 Weibo He.
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

#include "MatrixTrig.cuh"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    class device_obj<GEMV<M, V>> : public device_obj<RValueVector<GEMV<M, V>>> {
        using host_obj = GEMV<M, V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using RefM = add_device_obj<M>::type;
        using RefV = add_device_obj<V>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> inv;
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> rhs;
    public:
        device_obj(RefM inv, RefV rhs);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Vector auto& target) const;
        void assign_cublas(Vector auto& target) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return getRHS().getLength(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    device_obj<GEMV<M, V>>::device_obj(RefM inv, RefV rhs) : inv(asStruct(inv)), rhs(asStruct(rhs)) {}

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    void device_obj<GEMV<M, V>>::assign(Vector auto& target) const {
        assign_cublas(target);
    }

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    void device_obj<GEMV<M, V>>::assign_cublas(Vector auto& target) const {
        using Tm = decltype(std::declval<T>().toCUDA());
        getRHS().assign(target);

        const auto& mat = getLHS().getExpr().getExpr();
        constexpr auto TransA = MatrixMajor::isRowMatrix<decltype(mat)>() ? CUBLAS_OP_T : CUBLAS_OP_N;
        constexpr auto UploNoTrans = Traits<M>::Upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
        constexpr auto Uplo = []() consteval noexcept {
            if constexpr (TransA == CUBLAS_OP_T)
                return UploNoTrans == CUBLAS_FILL_MODE_UPPER ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
            return UploNoTrans;
        }();
        constexpr auto Diag = Traits<M>::Unit ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;

        const size_t n = getLength();
        const auto* a = reinterpret_cast<const Tm*>(mat.data_handle());
        const size_t lda = mat.getMajorStride();
        auto* x = reinterpret_cast<Tm*>(target.data_handle());

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (Base::isComplex()) {
            if constexpr (T::Prec == Float32)
                check(cublasCtrsv_64(ctx, Uplo, TransA, Diag, n, a, lda, x, 1));
            else
                check(cublasZtrsv_64(ctx, Uplo, TransA, Diag, n, a, lda, x, 1));
        }
        else {
            if constexpr (T::Prec == Float32)
                check(cublasStrsv_64(ctx, Uplo, TransA, Diag, n, a, lda, x, 1));
            else
                check(cublasDtrsv_64(ctx, Uplo, TransA, Diag, n, a, lda, x, 1));
        }
    }

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    auto&& device_obj<GEMV<M, V>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RefM>(self.inv.getDerived());
    }

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    auto&& device_obj<GEMV<M, V>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RefV>(self.rhs.getDerived());
    }
}
