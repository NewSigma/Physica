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

#include "../DiffDenseMatrix.cuh"

namespace Physica {
#define tparams Scalar T, int Order, int Option
#define DenseMatrix DenseMatrix<Diff<T, DiffMode::Reverse, Order>, Option>

    template<tparams>
    device_obj<DenseMatrix>::device_obj(size_t row, size_t col) : v(row, col), g(row, col) {}

    template<tparams>
    device_obj<DenseMatrix>::device_obj(size_t row, size_t col, T init) : v(row, col, init), g(row, col, 0) {}

    template<tparams>
    device_obj<DenseMatrix>::device_obj(ValueMatrix v_, GradMatrix g_) : v(std::move(v_)), g(std::move(g_)) {}

    template<tparams>
    template<Vector V>
    device_obj<DenseMatrix>::device_obj(const V& vec) requires(!ReverseDiff<V>)  : device_obj(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.assign(col);
    }

    template<tparams>
    template<Matrix M>
    device_obj<DenseMatrix>::device_obj(const M& mat) requires(!ReverseDiff<M>) : device_obj(mat.getRow(), mat.getCol()) {
        mat.assign(*this);
    }

    template<tparams>
    template<RNG R>
    inline void device_obj<DenseMatrix>::random_uniform() {
        *this = random_uniform<R>(getRow(), getCol());
    }

    template<tparams>
    template<RNG R>
    inline void device_obj<DenseMatrix>::random_normal() {
        *this = random_normal<R>(getRow(), getCol());
    }

    template<tparams>
    template<RNG R, class Distribution>
    inline void device_obj<DenseMatrix>::random_any(Distribution& dist) {
        *this = random_any<R, Distribution>(getRow(), getCol(), dist);
    }

    template<tparams>
    void device_obj<DenseMatrix>::resize(size_t row, size_t col) {
        v.resize(row, col);
        g.resize(row, col);
    }

    template<tparams>
    void device_obj<DenseMatrix>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<tparams>
    __host__ __device__ inline auto device_obj<DenseMatrix>::data_ptr(size_t row, size_t col) noexcept -> PtrTy {
        assert(row < getRow() && "[Error]: Index out of range");
        assert(col < getCol() && "[Error]: Index out of range");
        return PtrTy(v.data_ptr(row, col), g.data_ptr(row, col));
    }

    template<tparams>
    __host__ __device__ inline auto device_obj<DenseMatrix>::data_ptr(size_t row, size_t col) const noexcept -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<tparams>
    template<RNG R>
    inline auto device_obj<DenseMatrix>::random_uniform(size_t row, size_t col) -> This {
        return This(ValueMatrix::template random_uniform<R>(row, col));
    }

    template<tparams>
    template<RNG R>
    inline auto device_obj<DenseMatrix>::random_normal(size_t row, size_t col) -> This {
        return This(ValueMatrix::template random_normal<R>(row, col));
    }

    template<tparams>
    template<RNG R, class Distribution>
    inline auto device_obj<DenseMatrix>::random_any(size_t row, size_t col, Distribution& dist) -> This {
        return This(ValueMatrix::template random_any<R, Distribution>(row, col, dist));
    }

#undef tparams
#undef DenseMatrix
}
