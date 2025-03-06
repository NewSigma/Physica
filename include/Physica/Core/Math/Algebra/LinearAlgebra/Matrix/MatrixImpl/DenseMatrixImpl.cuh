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

#include "../DenseMatrix.cuh"

namespace Physica {
#define tparams Scalar T, int Option, size_t Row, size_t Col, class Allocator
#define DenseMatrix DenseMatrix<T, Option, Row, Col, Allocator>

    template<tparams>
    device_obj<DenseMatrix>::device_obj(const host_obj& mat) : device_obj(mat.getRow(), mat.getCol()) {
        mat.toDevice(*this);
    }

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix>::device_obj(size_t row, size_t col) : Storage(row, col) {}

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix>::device_obj(size_t row, size_t col, T value) : Storage(row, col, std::move(value)) {}

    template<tparams>
    template<Matrix M>
    device_obj<DenseMatrix>::device_obj(const device_obj<M>& mat) : device_obj(mat.getRow(), mat.getCol()) {
        mat.getDerived().assign(*this);
    }

    template<tparams>
    auto device_obj<DenseMatrix>::toHost() const -> host_obj {
        return host_obj(Storage::toHost());
    }

    template<tparams>
    auto device_obj<DenseMatrix>::toHostAsync() const-> host_obj {
        return host_obj(Storage::toHostAsync());
    }

    template<tparams>
    void device_obj<DenseMatrix>::toHost(host_obj& obj) const {
        Storage::toHost(obj);
    }

    template<tparams>
    void device_obj<DenseMatrix>::toHostAsync(host_obj& obj) const {
        Storage::toHostAsync(obj);
    }

    template<tparams>
    void device_obj<DenseMatrix>::zeros() {
        if constexpr (MatrixOption::isElementMatrix<This>())
            check(cudaMemsetAsync(data(), 0, getRow() * getCol() * sizeof(T), CUDAContext::getInstance()));
        else
            *this = Tv(0);
    }

    template<tparams>
    inline auto device_obj<DenseMatrix>::unitMatrix(size_t order) -> This {
        return host_obj::unitMatrix(order).toDevice();
    }

    template<tparams>
    template<RNG R>
    inline auto device_obj<DenseMatrix>::random_uniform(size_t row, size_t col) -> This {
        return host_obj::template random_uniform<R>(row, col).toDevice();
    }

    template<tparams>
    template<RNG R>
    inline auto device_obj<DenseMatrix>::random_normal(size_t row, size_t col) -> This {
        return host_obj::template random_normal<R>(row, col).toDevice();
    }

    template<tparams>
    template<RNG R, class Distribution>
    inline auto device_obj<DenseMatrix>::random_any(size_t row, size_t col, Distribution& dist) -> This {
        return host_obj::template random_any<R, Distribution>(row, col, dist).toDevice();
    }

#undef DenseMatrix
#undef tparams
}
