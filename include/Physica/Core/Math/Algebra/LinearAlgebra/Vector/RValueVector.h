/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixImpl/DenseMatrixOption.h"
#include "BestPacket.h"
#include "RVectorBlock.h"

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
        using ColMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector, SizeAtCompile, 1, MaxSizeAtCompile, 1>;
        using RowMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector, 1, SizeAtCompile, 1, MaxSizeAtCompile>;
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
        [[nodiscard]] PacketType packet(size_t index) const;
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const;
        [[nodiscard]] TransposeVector<Derived> transpose() const noexcept { return TransposeVector<Derived>(*this); }
        [[nodiscard]] ConjugateVector<Derived> conjugate() const noexcept { return ConjugateVector<Derived>(*this); }
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        [[nodiscard]] RealType norm() const;
        [[nodiscard]] RealType squaredNorm() const;
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
}

#include "RValueVectorImpl.h"
