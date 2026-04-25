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

#include "../RValueMatrixImpl.cuh"
#include "MatrixTrig.h"

namespace Physica {
    template<Matrix M, bool Upper, bool Unit>
    class device_obj<MatrixTrig<M, Upper, Unit>> : public device_obj<RValueMatrix<MatrixTrig<M, Upper, Unit>>> {
        using host_obj = MatrixTrig<M, Upper, Unit>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
        using typename Base::Trv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ device_obj(Ref mat);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;

        [[nodiscard]] auto lnAbsDet() const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getExpr().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getExpr().getCol(); }
    };

    template<Matrix M, bool Upper, bool Unit>
    __host__ __device__ device_obj<MatrixTrig<M, Upper, Unit>>::device_obj(Ref mat) : mat(asStruct(mat)) {}

    template<Matrix M, bool Upper, bool Unit>
    __device__ auto device_obj<MatrixTrig<M, Upper, Unit>>::calc(size_t row, size_t col) const -> T {
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
        return getExpr().calc(row, col);
    }

    template<Matrix M, bool Upper, bool Unit>
    auto device_obj<MatrixTrig<M, Upper, Unit>>::lnAbsDet() const {
        if constexpr (Unit)
            return Trv(0);
        else {
            assert(!Base::diag().prod().isZero() && "[Error]: Singular matrix");
            return ln(abs(Base::diag())).sum();
        }
    }

    template<Matrix M, bool Upper, bool Unit>
    auto device_obj<MatrixTrig<M, Upper, Unit>>::values(this auto&& self) noexcept {
        using Value = decltype(std::forward<decltype(self)>(self).getExpr().values())::host_obj;
        return device_obj<MatrixTrig<Value, Upper, Unit>>(std::forward<decltype(self)>(self).getExpr().values());
    }

    template<Matrix M, bool Upper, bool Unit>
    __host__ __device__ auto&& device_obj<MatrixTrig<M, Upper, Unit>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }
}

namespace Physica {
    template<Matrix M, bool Upper, bool Unit>
    class Traits<device_obj<MatrixTrig<M, Upper, Unit>>> : public Traits<MatrixTrig<M, Upper, Unit>> {};
}

#include "GEMM.cuh"
#include "Inverse.cuh"
#include "InvGEMM.cuh"
