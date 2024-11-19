/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh>
#include "MatrixImpl/ContinuousMatrix.cuh"
#include "DenseMatrixImpl/DenseMatrixStorage.cuh"

namespace Physica::Core {
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    class device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>
            : public device_obj<ContinuousMatrix<DenseMatrix<T, Option, Row, Col, Allocator>>>
            , public device_obj<DenseMatrixStorage<T, Option, Row, Col, Allocator>> {
        using host_obj = DenseMatrix<T, Option, Row, Col, Allocator>;
        using host_storage = typename host_obj::Storage;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousMatrix<DenseMatrix<T, Option, Row, Col, Allocator>>>;
        using Storage = device_obj<DenseMatrixStorage<T, Option, Row, Col, Allocator>>;
    public:
        device_obj() = default;
        __host__ __device__ device_obj(size_t row, size_t col);
        __host__ __device__ device_obj(size_t row, size_t col, T value);
        device_obj(const host_obj& mat);
        template<Matrix M>
        device_obj(const device_obj<M>& mat);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Storage::operator();
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Storage::toHost()); }
        using Storage::resize;
        using Storage::swap;
        /* Getters */
        using Base::data;
        using Storage::data_ptr;
        using Storage::getRow;
        using Storage::getCol;
        /* Static members */
        [[nodiscard]] static inline device_obj unitMatrix(size_t order);
        template<class RandomType>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t col, RandomType& gen);
        template<class RandomType>
        [[nodiscard]] inline static This random_normal(size_t row, size_t col, RandomType& gen);
        template<class Distribution, class RandomType>
        [[nodiscard]] inline static This random_any(size_t row, size_t col, Distribution& dist, RandomType& gen);
    };

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(
            size_t row, size_t col) : Storage(row, col) {}

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(
            size_t row, size_t col, T value) : Storage(row, col, std::move(value)) {}

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(const host_obj& mat)
            : Storage(static_cast<const host_storage&>(mat).toDevice()) {}

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<Matrix M>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::device_obj(
            const device_obj<M>& mat) : device_obj(mat.getRow(), mat.getCol()) {
        mat.getDerived().assignTo(*this);
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    inline device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::unitMatrix(size_t order) {
        return host_obj::unitMatrix(order).toDevice();
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<class RandomType>
    inline device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_uniform(
            size_t row, size_t col, RandomType& gen) {
        return host_obj::random_uniform(row, col, gen).toDevice();
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<class RandomType>
    inline device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_normal(
            size_t row, size_t col, RandomType& gen) {
        return host_obj::random_normal(row, col, gen).toDevice();
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    template<class Distribution, class RandomType>
    inline device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>::random_any(
            size_t row, size_t col, Distribution& dist, RandomType& gen) {
        return host_obj::random_any(row, col, dist, gen).toDevice();
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    inline auto DenseMatrix<T, Option, Row, Col, Allocator>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    inline auto DenseMatrix<T, Option, Row, Col, Allocator>::toDeviceAsync() const {
        device_obj<This> result(getRow(), getCol());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Option, Row, Col, Allocator>::toDevice(device_obj<This>& obj) const {
        obj.resize(getRow(), getCol());
        Storage::toDevice(obj);
    }

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Option, Row, Col, Allocator>::toDeviceAsync(device_obj<This>& obj) const {
        obj.resize(getRow(), getCol());
        Storage::toDeviceAsync(obj);
    }
}

namespace Physica {
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    class Traits<Core::device_obj<DenseMatrix<T, Option, Row, Col, Allocator>>>
            : public Traits<DenseMatrix<T, Option, Row, Col, Allocator>> {};
}
