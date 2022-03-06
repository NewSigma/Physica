/*
 * Copyright 2020-2021 WeiBo He.
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
#include "ContinuousVector.h"
#include "VectorBlock.h"
#include "VectorExpression.h"
#include "BestPacket.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixImpl/DenseMatrixOption.h"

namespace Physica::Core {
    template<class T = MultiScalar, size_t Length = Dynamic, size_t MaxLength = Length>
    class Vector;

    namespace Internal {
        template<class T, size_t Length, size_t MaxLength>
        class Traits<Vector<T, Length, MaxLength>> {
        public:
            using ScalarType = T;
            constexpr static size_t SizeAtCompile = Length;
            constexpr static size_t MaxSizeAtCompile = MaxLength;
        };
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrix;
    /**
     * T must be either Scalar or ComplexScalar.
     */
    template<class T, size_t Length, size_t MaxLength>
    class Vector : public ContinuousVector<Vector<T, Length, MaxLength>>, public Utils::Array<T, Length, MaxLength> {
        static_assert(Length == Dynamic || Length == MaxLength, "MaxLength of fixed vector must equals to its length.");
        using Storage = Utils::Array<T, Length, MaxLength>;
    public:
        using Base = ContinuousVector<Vector<T, Length, MaxLength>>;
        using typename Base::ColMatrix;
        using typename Base::RowMatrix;
    public:
        using Storage::Storage;
        Vector() = default;
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
        Vector& toOpposite();
        void toUnit();
        /* Getters */
        using Storage::getLength;
        /* Helpers */
        using Storage::swap;
        static Vector Zeros(size_t len);
        static Vector randomVector(size_t len);
        static Vector randomVector(const Vector& v1, const Vector& v2);
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

#include "VectorImpl.h"
