/*
 * Copyright 2024 WeiBo He.
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

#include "DiffDenseMatrix.h"

namespace Physica::Core {
    namespace Internal {
        template<class T, int Option>
        class Traits<device_obj<Differentiable<DenseMatrix<T, Option>, DiffMode::Reverse>>> : public Traits<DenseMatrix<T, Option>> {
        public:
            using ScalarType = device_obj<Differentiable<T, DiffMode::Reverse>>;
        };
    }

    template<class PlainScalar, int Option>
    class device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>>
            : public device_obj<RValueMatrix<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>>>
            , public DenseMatrixDim<device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>>> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using host_obj = Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using SegmentType = device_obj<typename host_obj::SegmentType>;
        using Dim = DenseMatrixDim<This>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        PlainStruct<SegmentType> traceSeg;
    public:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        /* Operations */
        void swap(device_obj& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getRow;
        using Dim::getColumn;
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return traceSeg.getLength(); }
    private:
        [[nodiscard]] SegmentType& getTraceSegment() noexcept { return traceSeg.getDerived(); }
    };
}
