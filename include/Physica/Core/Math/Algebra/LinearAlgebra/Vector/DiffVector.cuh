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

#include "Physica/Core/MultiPrecision/Differentiable.cuh"
#include "DiffVector.h"
#include "Vector.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class T>
        class Traits<device_obj<Differentiable<Vector<T>, DiffMode::Reverse>>> : public Traits<Vector<T>> {
        public:
            using ScalarType = device_obj<Differentiable<T, DiffMode::Reverse>>;
        };
    }

    template<class PlainScalar>
    class device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>
            : public device_obj<RValueVector<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using VectorType = Vector<PlainScalar>;
        using host_obj = Differentiable<VectorType, DiffMode::Reverse>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using TracerType = device_obj<typename host_obj::TracerType>;
        using SegmentType = device_obj<typename host_obj::SegmentType>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        PlainStruct<SegmentType> traceSeg;
    public:
        device_obj() = default;
        device_obj(const VectorType& values);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        /* Operations */
        void swap(device_obj& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getTraceSegment().getLength(); }
    private:
        [[nodiscard]] SegmentType& getTraceSegment() noexcept { return traceSeg.getDerived(); }
    }

    template<class PlainScalar>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::device_obj(const VectorType& values)
            : traceSeg(asStruct(TracerType::getInstance().pushSegment(values))) {}

    template<class PlainScalar>
    __host__ __device__ inline typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        return getTraceSegment()[index];
    }
}
