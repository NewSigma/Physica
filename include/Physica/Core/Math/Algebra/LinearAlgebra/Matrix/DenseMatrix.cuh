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

#include "MatrixImpl/ContinuousMatrix.cuh"
#include "DenseMatrixImpl/DenseMatrixStorage.cuh"
#include "MatrixProduct/MatrixProduct.cuh"

namespace Physica::Core {
    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
            : public device_obj<ContinuousMatrix<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>>
            , public device_obj<DenseMatrixStorage<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>> {
        using host_obj = DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>;
        using host_storage = typename host_obj::Storage;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using Storage = device_obj<DenseMatrixStorage<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>;
    public:
        device_obj() = default;
        __host__ __device__ device_obj(size_t row, size_t column);
        __host__ __device__ device_obj(size_t row, size_t column, T value);
        device_obj(const host_obj& mat);
        template<class OtherMatrix>
        device_obj(const device_obj<RValueMatrix<OtherMatrix>>& mat);
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
        using Storage::getRow;
        using Storage::getColumn;
        using Storage::data_ptr;
        /* Static members */
        [[nodiscard]] static inline device_obj unitMatrix(size_t order);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t column, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t row, size_t column, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t row, size_t column, Distribution& dist, RandomGenerator& gen);
    };

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::device_obj(
            size_t row, size_t column) : Storage(row, column) {}

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::device_obj(
            size_t row, size_t column, T value) : Storage(row, column, std::move(value)) {}

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::device_obj(const host_obj& mat)
            : Storage(static_cast<const host_storage&>(mat).toDevice()) {}

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class OtherMatrix>
    device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::device_obj(
            const device_obj<RValueMatrix<OtherMatrix>>& mat) : device_obj(mat.getRow(), mat.getColumn()) {
        mat.getDerived().assignTo(*this);
    }

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::unitMatrix(size_t order) {
        return host_obj::unitMatrix(order).toDevice();
    }

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class RandomGenerator>
    inline device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::random_uniform(
            size_t row, size_t column, RandomGenerator& gen) {
        return host_obj::random_uniform(row, column, gen).toDevice();
    }

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class RandomGenerator>
    inline device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::random_normal(
            size_t row, size_t column, RandomGenerator& gen) {
        return host_obj::random_normal(row, column, gen).toDevice();
    }

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class Distribution, class RandomGenerator>
    inline device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
    device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>::random_any(
            size_t row, size_t column, Distribution& dist, RandomGenerator& gen) {
        return host_obj::random_any(row, column, dist, gen).toDevice();
    }

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
    DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>::toDevice() const {
        return device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>(*this);
    }

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>::toDevice(device_obj<This>& obj) const {
        obj.resize(getRow(), getColumn());
        Storage::toDevice(obj);
    }
}

namespace Physica {
    using namespace Core;

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class Traits<Core::device_obj<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>>
            : public Traits<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>> {};
}
