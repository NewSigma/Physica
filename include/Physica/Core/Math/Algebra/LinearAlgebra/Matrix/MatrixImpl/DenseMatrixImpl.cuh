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

    template<tparams>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(const host_obj& mat)
            : device_obj(mat.getRow(), mat.getCol()) {
        mat.toDevice(*this);
    }

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(size_t row, size_t col)
            : Storage(row, col) {}

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(size_t row, size_t col, T value)
            : Storage(row, col) {
        *this = value;
    }

    template<tparams>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(const Matrix auto& mat) : device_obj(mat.getRow(), mat.getCol()) {
        mat.assign(*this);
    }

    template<tparams>
    __host__ __device__ void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::resize(const Matrix auto& m) {
        Storage::resize(m.getRow(), m.getCol());
    }

    template<tparams>
    auto device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::toHost() const -> host_obj {
        return host_obj(Storage::toHost());
    }

    template<tparams>
    auto device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::toHostAsync() const -> host_obj {
        return host_obj(Storage::toHostAsync());
    }

    template<tparams>
    void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::toHost(host_obj& obj) const {
        Storage::toHost(obj);
    }

    template<tparams>
    void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::toHostAsync(host_obj& obj) const {
        Storage::toHostAsync(obj);
    }

    template<tparams>
    void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::zeros() {
        check(cudaMemsetAsync(data(), 0, getRow() * getCol() * sizeof(T), CUDAContext::getInstance()));
    }

    template<tparams>
    template<RNG R>
    void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_uniform() {
        check(curandSetStream(R::getInstance(), CUDAContext::getInstance()));
        if constexpr (T::Prec == Float32)
            check(curandGenerateUniform(R::getInstance(), (Tm*)data(), getRow() * getCol()));
        else if constexpr (T::Prec == Float64)
            check(curandGenerateUniformDouble(R::getInstance(), (Tm*)data(), getRow() * getCol()));
        else
            host_obj::template random_uniform<R>(getRow(), getCol()).toDeviceAsync(*this);
    }

    template<tparams>
    template<RNG R>
    void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_normal() {
        check(curandSetStream(R::getInstance(), CUDAContext::getInstance()));
        if constexpr (T::Prec == Float32)
            check(curandGenerateNormal(R::getInstance(), (Tm*)data(), getRow() * getCol(), 0, 1));
        else if constexpr (T::Prec == Float64)
            check(curandGenerateNormalDouble(R::getInstance(), (Tm*)data(), getRow() * getCol(), 0, 1));
        else
            *this = device_obj<DenseMatrix<float32, Option, Row, Col, Allocator>>::template random_normal<R>(getRow(), getCol());
    }

    template<tparams>
    template<RNG R>
    void device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_any(auto& distribution) {
        host_obj::template random_any<R>(getRow(), getCol(), distribution).toDeviceAsync(*this);
    }

    template<tparams>
    auto device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::unitMatrix(size_t order) -> This {
        return host_obj::unitMatrix(order).toDevice();
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_uniform(size_t row, size_t col) -> This {
        This result(row, col);
        result.template random_uniform<R>();
        return result;
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_normal(size_t row, size_t col) -> This {
        This result(row, col);
        result.template random_normal<R>();
        return result;
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_any(size_t row, size_t col, auto& distribution) -> This {
        This result(row, col);
        result.template random_any<R>(distribution);
        return result;
    }

#undef tparams
}
