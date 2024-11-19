/*
 * Copyright 2020-2024 Weibo He.
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

#include <Physica/Core/MultiPrecision/Real.h>
#include <Physica/Core/Utils/Container/Array.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h>
#include "VectorImpl/ContinuousVector.h"

namespace Physica::Core {
    template<Scalar T, size_t Length = Dynamic, class Allocator = HostAllocator<T, alignof(typename BestPacket<T, Length>::Type)>>
    class DenseVector : public ContinuousVector<DenseVector<T, Length, Allocator>>, public Array<T, Length, Allocator> {
        constexpr static size_t DefaultAlign = alignof(typename BestPacket<T, Length>::Type);
        static_assert(std::allocator_traits<Allocator>::Align % DefaultAlign == 0, "[Error]: Bad alignment for SIMD");
        using This = DenseVector<T, Length, Allocator>;
    public:
        using Base = ContinuousVector<This>;
        using Storage = Array<T, Length, Allocator>;
        using device_obj_type = device_obj<This>;
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using typename Base::ColMatrix;
        using typename Base::RowMatrix;
        using Base::SizeAtCompile;
        using Base::isReverseDiff;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    public:
        DenseVector() = default;
        explicit DenseVector(size_t length);
        DenseVector(size_t length, const T& init);
        DenseVector(std::initializer_list<T> list);
        explicit DenseVector(Storage array) noexcept;
        template<Vector V>
        DenseVector(const V& v);
        DenseVector(const This&) = default;
        DenseVector(This&&) noexcept = default;
        ~DenseVector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        using Storage::resize;
        using Base::random_uniform;
        using Base::random_normal;
        using Base::random_any;
        [[nodiscard]] This copy() const;
        [[nodiscard]] inline device_obj<This> toDevice() const;
        using Base::toDevice;
        using Base::toDeviceAsync;
        using Storage::swap;
        /* Getters */
        using Storage::getLength;
        using Storage::data;
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) { return data() + index; }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const { return data() + index; }
        /* Static members */
        [[nodiscard]] static This zeros(size_t len);
        template<class RandomGenerator>
        [[nodiscard]] static This random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static This random_uniform(const This& v1, const This& v2, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static This random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] static This random_any(size_t len, Distribution& dist, RandomGenerator& gen);
        [[nodiscard]] static This linspace(T from, T to, size_t count);
    };

    template<Scalar T> using Vector1D = DenseVector<T, 1>;
    template<Scalar T> using Vector2D = DenseVector<T, 2>;
    template<Scalar T> using Vector3D = DenseVector<T, 3>;
    template<Scalar T> using Vector4D = DenseVector<T, 4>;
    template<Scalar T> using VectorND = DenseVector<T>;
}

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    class Traits<Core::DenseVector<T, Length, Allocator>> {
        static_assert(!T::isForwardDiff, "[Error]: Use diffable vector instead");
    public:
        using ScalarType = T;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}

namespace std {
    template<Physica::Core::Scalar T, size_t Length>
    inline void swap(Physica::Core::DenseVector<T, Length>& __restrict v1, Physica::Core::DenseVector<T, Length>& __restrict v2) noexcept {
        v1.swap(v2);
    }
}

#include "VectorImpl/DenseVectorImpl.h"
