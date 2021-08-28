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
#include "LValueVector.h"
#include "VectorBlock.h"
#include "VectorExpression.h"
#include "CrossProduct.h"
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
    class Vector : public LValueVector<Vector<T, Length, MaxLength>>, public Utils::Array<T, Length, MaxLength> {
        static_assert(Length == Dynamic || Length == MaxLength, "MaxLength of fixed vector must equals to its length.");
        using Storage = Utils::Array<T, Length, MaxLength>;
    public:
        using Base = LValueVector<Vector<T, Length, MaxLength>>;
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
        Vector& operator=(const Vector&) = default;
        Vector& operator=(Vector&&) noexcept = default;
        using Storage::operator[];
        /* Operations */
        Vector& toOpposite();
        void toUnit();
        template<class OtherVector>
        [[nodiscard]] inline CrossProduct<Vector, OtherVector> crossProduct(const RValueVector<OtherVector>& v) const noexcept;
        ColMatrix moveToColMatrix();
        RowMatrix moveToRowMatrix();
        /* Getters */
        using Storage::getLength;
        [[nodiscard]] bool isZero() const;
        [[nodiscard]] T max() const;
        [[nodiscard]] T min() const;
        /* Helpers */
        static Vector<T> Zeros(size_t len);
        static Vector<T> randomVector(size_t len);
        static Vector simplyMultiply(const Vector& v1, const Vector& v2);
    private:
        template<class Derived>
        friend class Internal::VectorExpressionHelper;
    };
    /* Operators */
    template<class T, size_t Length, size_t MaxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, Length, MaxLength>& v);
}

#include "VectorImpl.h"
