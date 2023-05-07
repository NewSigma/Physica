/*
 * Copyright 2021-2023 WeiBo He.
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

#include "Physica/Utils/Template/CRTPBase.h"
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h"
#include "VectorImpl/BestPacket.h"
#include "VectorImpl/RVectorBlock.h"

namespace Physica::Core {
    template<class VectorType> class TransposeVector;
    template<class VectorType> class ConjugateVector;
    template<class AnyVector1, class AnyVector2> class CrossProduct;
    template<class VectorType> class FormatedVector;

    namespace Internal {
        template<class T> class Traits;
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn> class DenseMatrix;

    template<class Derived> class LValueVector;
    /**
     * \class RValueVector is base class of vectors that can be assigned to \class LValueVector
     * but other vectors cannot be assigned to this class.
     * In other words, you cannot take the address of elements in the vector.
     */
    template<class Derived>
    class RValueVector : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        constexpr static size_t SizeAtCompile = Internal::Traits<Derived>::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Internal::Traits<Derived>::MaxSizeAtCompile;
        using ColMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector, SizeAtCompile, 1, MaxSizeAtCompile, 1>;
        using RowMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector, 1, SizeAtCompile, 1, MaxSizeAtCompile>;
        constexpr static bool isComplex = ScalarType::isComplex;
    private:
        using RealType = typename ScalarType::RealType;
    public:
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueVector<OtherDerived>& v) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] FormatedVector<Derived> format() const;
        template<class PacketType>
        [[nodiscard]] inline PacketType packet(size_t index) const;
        template<class PacketType>
        [[nodiscard]] inline PacketType packetPartial(size_t index) const;
        [[nodiscard]] TransposeVector<Derived> transpose() const noexcept { return TransposeVector<Derived>(*this); }
        [[nodiscard]] ConjugateVector<Derived> conjugate() const noexcept { return ConjugateVector<Derived>(*this); }
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        [[nodiscard]] inline RealType norm() const;
        [[nodiscard]] inline RealType squaredNorm() const;
        [[nodiscard]] ScalarType max() const;
        [[nodiscard]] ScalarType min() const;
        [[nodiscard]] ScalarType sum() const;
        template<class OtherDerived>
        [[nodiscard]] inline CrossProduct<Derived, OtherDerived> crossProduct(const RValueVector<OtherDerived>& v) const noexcept;
        template<class OtherDerived>
        [[nodiscard]] ScalarType angleTo(const RValueVector<OtherDerived>& v) const noexcept;
        RVectorBlock<Derived> head(size_t to) { return RVectorBlock<Derived>(Base::getDerived(), 0, to); }
        const RVectorBlock<Derived> head(size_t to) const { return RVectorBlock<Derived>(Base::getConstCastDerived(), 0, to); }
        RVectorBlock<Derived> tail(size_t from) { return RVectorBlock<Derived>(Base::getDerived(), from); }
        const RVectorBlock<Derived> tail(size_t from) const { return RVectorBlock<Derived>(Base::getConstCastDerived(), from); }
        RVectorBlock<Derived> segment(size_t from, size_t to) { return RVectorBlock<Derived>(Base::getDerived(), from, to); }
        const RVectorBlock<Derived> segment(size_t from, size_t to) const { return RVectorBlock<Derived>(Base::getConstCastDerived(), from, to); }
    };

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2);

    template<class VectorType1, class VectorType2>
    bool vectorNear(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2, double precision);

    template<class VectorType>
    std::ostream& operator<<(std::ostream& os, const RValueVector<VectorType>& v);
}

#include "VectorImpl/RValueVectorImpl.h"
#include "VectorImpl/VectorExpression.h"
