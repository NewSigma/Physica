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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<class M, bool ReduceCol>
    class device_obj<MatrixSum<M, ReduceCol>> : public device_obj<RValueVector<MatrixSum<M, ReduceCol>>> {
        using host_obj = MatrixSum<M, ReduceCol>;
        using This = device_obj<MatrixSum<M, ReduceCol>>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<M>> mat;
    public:
        __host__ __device__ device_obj(const device_obj<M>& mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index) const;

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ const device_obj<M>& getExpr() const noexcept { return mat.getDerived(); }
    };

    template<class M, bool ReduceCol>
    __device__ auto device_obj<MatrixSum<M, ReduceCol>>::calc(size_t index) const -> T {
        if constexpr (ReduceCol)
            return mat.getDerived().row(index).sum();
        else
            return mat.getDerived().col(index).sum();
    }

    template<class M, bool ReduceCol>
    __device__ auto device_obj<MatrixSum<M, ReduceCol>>::calc_value(size_t index) const -> Tv {
        if constexpr (ReduceCol)
            return mat.getDerived().values().row(index).sum();
        else
            return mat.getDerived().values().col(index).sum();
    }

    template<class M, bool ReduceCol>
    void device_obj<MatrixSum<M, ReduceCol>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        if constexpr (ReduceCol)
            mat.getDerived().reverse(grad);
        else
            mat.getDerived().transpose().reverse(grad);
    }

    template<class M, bool ReduceCol>
    __host__ __device__ size_t device_obj<MatrixSum<M, ReduceCol>>::getLength() const noexcept {
        if constexpr (ReduceCol)
            return mat.getDerived().getRow();
        else
            return mat.getDerived().getCol();
    }
}

namespace Physica {
    template<Matrix M, bool ReduceCol>
    class Traits<device_obj<MatrixSum<M, ReduceCol>>> : public Traits<MatrixSum<M, ReduceCol>> {};
}
