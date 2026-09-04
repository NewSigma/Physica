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

#include "../DenseMatrix.cuh"

namespace Physica {
#define tparams Scalar T, int Major, size_t Row, size_t Col, class Allocator

    template<tparams>
    device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(Storage storage) noexcept
            : storage(std::move(storage)) {}

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(size_t order)
            : storage(order) {}

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(size_t row, size_t col)
            : storage(row, col) {}

    template<tparams>
    __host__ __device__ device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(size_t row, size_t col, T value)
            : storage(row, col) {
        *this = value;
    }

    template<tparams>
    device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(const Matrix auto& mat) : device_obj(mat.getRow(), mat.getCol()) {
        mat.assign(*this);
    }

    template<tparams>
    device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(const Vector auto& vec) : device_obj(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.assign(col);
    }

    template<tparams>
    device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::device_obj(const host_obj& mat)
            : device_obj(mat.getRow(), mat.getCol()) {
        mat.toDevice(*this);
    }

    template<tparams>
    __host__ __device__ void device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::resize(this auto& self, size_t row, size_t col, auto&&... args) {
        self.storage.resize(row, col, std::forward<decltype(args)>(args)...);
    }

    template<tparams>
    void device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::reserve(this auto& self, size_t size) {
        self.storage.reserve(size);
    }

    template<tparams>
    auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::toHost() const -> host_obj {
        return host_obj(storage.toHost());
    }

    template<tparams>
    auto&& device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::flatten(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.storage).asArray();
    }

    template<tparams>
    auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::toHostAsync() const -> host_obj {
        return host_obj(storage.toHostAsync());
    }

    template<tparams>
    void device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::zeros(this auto& self) noexcept {
        self.storage.zeros();
    }

    template<tparams>
    void device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        storage.swap(obj.storage);
    }

    template<tparams>
    __host__ __device__ auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::data(this auto&& self) noexcept {
        return self.storage.data();
    }

    template<tparams>
    __host__ __device__ auto&& device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::asArray(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.storage).asArray();
    }

    template<tparams>
    auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::identity(size_t order) -> This {
        return host_obj::identity(order).toDevice();
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::random_uniform(size_t row, size_t col) -> This {
        This result(row, col);
        result.template random_uniform<R>();
        return result;
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::random_normal(size_t row, size_t col) -> This {
        This result(row, col);
        result.template random_normal<R>();
        return result;
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>::random_any(size_t row, size_t col, auto& distribution) -> This {
        This result(row, col);
        result.template random_any<R>(distribution);
        return result;
    }

#undef tparams
}
