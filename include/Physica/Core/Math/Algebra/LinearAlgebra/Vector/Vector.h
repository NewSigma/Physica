/*
 * Copyright 2020-2023 WeiBo He.
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

#include <iosfwd>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Container/Array/Array.h"
#include "VectorImpl/ContinuousVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h"

namespace Physica::Core {
    template<class T = MultiScalar, size_t Length = Dynamic, size_t MaxLength = Length, class Allocator = Utils::HostAllocator<T>>
    class Vector;

    namespace Internal {
        template<class T, size_t Length, size_t MaxLength, class Allocator>
        class Traits<Vector<T, Length, MaxLength, Allocator>> {
        public:
            using ScalarType = T;
            constexpr static size_t SizeAtCompile = Length;
            constexpr static size_t MaxSizeAtCompile = MaxLength;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class T, size_t Length, size_t MaxLength, class Allocator>
    class Vector : public ContinuousVector<Vector<T, Length, MaxLength, Allocator>>, public Utils::Array<T, Length, MaxLength, Allocator> {
        static_assert(Length == Dynamic || Length == MaxLength, "MaxLength of fixed vector must equals to its length.");
        using This = Vector<T, Length, MaxLength, Allocator>;
        using Storage = Utils::Array<T, Length, MaxLength, Allocator>;
    public:
        using Base = ContinuousVector<This>;
        using device_obj_type = device_obj<This>;
        using typename Base::ColMatrix;
        using typename Base::RowMatrix;
    public:
        using Storage::Storage;
        Vector() = default;
        Vector(Storage array) noexcept;
        template<class Derived>
        Vector(const RValueVector<Derived>& v);
        Vector(const Vector&) = default;
        Vector(Vector&&) noexcept = default;
        ~Vector() = default;
        /* Operators */
        Vector& operator=(Vector v) noexcept;
        using Base::operator=;
        using Storage::operator[];
        /* Operations */
        using Storage::resize;
        using Base::random_uniform;
        using Base::random_normal;
        Vector& toOpposite();
        void toUnit();
        [[nodiscard]] inline device_obj<This> toDevice() const;
        inline void toDevice(device_obj<This>& obj) const;
        /* Getters */
        using Storage::getLength;
        using Storage::data;
        [[nodiscard]] Vector reverse() const;
        /* Helpers */
        using Storage::swap;
        static Vector Zeros(size_t len);
        template<class RandomGenerator>
        static Vector random(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        static Vector random(const Vector& v1, const Vector& v2, RandomGenerator& gen);
        template<class RandomGenerator>
        static Vector random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        static Vector random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        static Vector random_any(size_t len, Distribution& dist, RandomGenerator& gen);
        static Vector linspace(T from, T to, size_t count);
    };

    namespace Internal {
        //TODO: extend RValueVector for optimized performance
        template<class VectorType>
        class RealVectorReturnType {
            using ComplexType = typename VectorType::ScalarType;
            using ScalarType = Scalar<Traits<ComplexType>::option, Traits<ComplexType>::errorTrack>;
            constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = VectorType::MaxSizeAtCompile;
        public:
            using Type = Vector<ScalarType, SizeAtCompile, MaxSizeAtCompile>;
        };
    }

    template<class VectorType>
    [[nodiscard]] typename Internal::RealVectorReturnType<VectorType>::Type toRealVector(const RValueVector<VectorType>& v) {
        using ResultType = typename Internal::RealVectorReturnType<VectorType>::Type;
        ResultType result = ResultType(v.getLength());
        for (size_t i = 0; i < v.getLength(); ++i)
            result[i] = v.calc(i).getReal();
        return result;
    }

    template<class VectorType>
    [[nodiscard]] typename Internal::RealVectorReturnType<VectorType>::Type toImagVector(const RValueVector<VectorType>& v) {
        using ResultType = typename Internal::RealVectorReturnType<VectorType>::Type;
        ResultType result = ResultType(v.getLength());
        for (size_t i = 0; i < v.getLength(); ++i)
            result[i] = v.calc(i).getImag();
        return result;
    }

    template<class VectorType>
    [[nodiscard]] typename Internal::RealVectorReturnType<VectorType>::Type toNormVector(const RValueVector<VectorType>& v) {
        using ResultType = typename Internal::RealVectorReturnType<VectorType>::Type;
        ResultType result = ResultType(v.getLength());
        for (size_t i = 0; i < v.getLength(); ++i)
            result[i] = v.calc(i).norm();
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    inline void swap(Vector<T, Length, MaxLength>& v1, Vector<T, Length, MaxLength>& v2) noexcept {
        v1.swap(v2);
    }
}

#include "VectorImpl/VectorImpl.h"
