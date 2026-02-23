/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Utils/Container/Array2D.cuh"
#include "MatrixImpl/CompactMatrix.cuh"
#include "DenseMatrix.h"

namespace Physica {
    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    class device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>
            : public device_obj<CompactMatrix<DenseMatrix<T, Major, Row, Col, Allocator>>>
            , public CRCoro<device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>>
            , public device_obj<Array2D<T, Major, Row, Col, Allocator>> {
        static_assert(!Diffable<T>, "[Error]: Use diffable matrix instead");
        using host_obj = DenseMatrix<T, Major, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactMatrix<DenseMatrix<T, Major, Row, Col, Allocator>>>;
        using Coro = CRCoro<This>;
        using Storage = device_obj<Array2D<T, Major, Row, Col, Allocator>>;
    protected:
        using typename Base::Tv;
    public:
        device_obj() = default;
        explicit device_obj(const host_obj& mat);
        explicit __host__ __device__ device_obj(size_t order);
        __host__ __device__ device_obj(size_t row, size_t col);
        __host__ __device__ device_obj(size_t row, size_t col, T value);
        device_obj(const Matrix auto& mat);
        device_obj(const Vector auto& vec);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Storage::operator[];
        /* Operations */
        __host__ __device__ void resize(const Matrix auto& m);
        using Storage::resize;
        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        using Base::toHost;
        using Base::toHostAsync;

        using Base::transpose;
        [[nodiscard]] auto&& flatten(this auto&&) noexcept;

        using Base::zeros;
        using Base::random_uniform;
        using Base::random_normal;
        using Base::random_any;

        using Storage::swap;
        /* Getters */
        using Storage::data;
        using Storage::data_ptr;
        using Storage::getRow;
        using Storage::getCol;
        using Storage::getSize;
        /* Static members */
        [[nodiscard]] static device_obj identity(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t row, size_t col, auto& distribution);
    };

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::toDeviceAsync() const {
        device_obj<This> result(getRow(), getCol());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }
}

namespace Physica {
    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    class Traits<device_obj<DenseMatrix<T, Major, Row, Col, Allocator>>>
            : public Traits<DenseMatrix<T, Major, Row, Col, Allocator>> {};
}

#include "MatrixImpl/DenseMatrixImpl.cuh"
