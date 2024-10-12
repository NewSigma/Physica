/*
 * Copyright 2021-2024 Weibo He.
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

#include <Physica/CRTPBase.h>
#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Core/MultiPrecision/Complex.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h>
#include <Physica/Core/Parallel/Executor/SequentialExecutor.h>
#include "RValueVectorImpl/RVectorBlock.h"

namespace Physica::Core {
    template<class Derived> class LValueVector;
    template<class Derived> class ContinuousVector;
    template<class Derived> class ContinuousMatrix;
    template<class T, int Option, size_t Row, size_t Column, class Allocator> class DenseMatrix;
    template<class VectorType> class TransposeVector;
    template<class VectorType> class ConjugateVector;
    template<class VectorType> class HermiteVector;
    template<class VectorType, int MatrixMajor, size_t Row, size_t Column> class ReshapedVector;
    template<class VectorType> class FormatedVector;
    template<class VectorType> class ReverseVector;
    template<class AnyVector1, class AnyVector2> class CrossProduct;
    template<class Derived> class ContinuousVector;
    template<class Derived> class RValueMatrix;

    template<class T>
    struct is_vector {
        constexpr static bool value = std::is_base_of<RValueVector<T>, T>::value;
    };

    template<class T>
    struct is_continuous {
        constexpr static bool value = std::is_base_of<ContinuousVector<T>, T>::value || std::is_base_of<ContinuousMatrix<T>, T>::value;
    };

    namespace Internal {
        template<class VectorType1, class VectorType2 = VectorType1>
        class EnableSIMD {
            using ScalarType = typename VectorType1::ScalarType;
            constexpr static bool isSameScalar = std::is_same<ScalarType, typename VectorType2::ScalarType>::value;
        public:
            constexpr static bool value = isSameScalar && BestPacket<ScalarType, VectorType1::SizeAtCompile>::Size > 1;
        };
    }

    template<class VectorType1, class VectorType2 = VectorType1>
    struct EnableMKL {
        constexpr static bool value = HasMKL() && is_continuous<VectorType1>::value && is_continuous<VectorType2>::value;
    };
    /**
     * \class RValueVector is a base class for vectors. You cannot take the address of elements in an RValueVector.
     * An RValueVector can be assigned to an LValueVector, but no other vector classes can be assigned to an RValueVector.
     */
    template<class Derived>
    class RValueVector : public CRTPBase<RValueVector<Derived>> {
        static_assert(!std::is_const<Derived>::value, "[Error]: A common mistake, const is unnecessary");
        static_assert(!std::is_volatile<Derived>::value, "[Error]: A common mistake, volatile is unnecessary");
        using This = RValueVector<Derived>;
        using Base = CRTPBase<This>;
        template<size_t Length>
        using BlockType = RVectorBlock<Derived, Length>;
    public:
        using ScalarType = typename Traits<Derived>::ScalarType;
        using PlainScalar = typename ScalarType::PlainScalar;
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        using PacketType = typename BestPacket<ScalarType, SizeAtCompile>::Type;
        using ColMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector, SizeAtCompile, 1, HostAllocator<ScalarType>>;
        using RowMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector, 1, SizeAtCompile, HostAllocator<ScalarType>>;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
    private:
        using RealType = typename ScalarType::RealType;
        constexpr static bool isExpression = !std::is_base_of<LValueVector<Derived>, Derived>::value;
    public:
        ~RValueVector() = default;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& v) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        template<class AnyPacket>
        [[nodiscard]] inline AnyPacket packet(size_t index) const;
        template<class AnyPacket>
        [[nodiscard]] inline AnyPacket packetPartial(size_t index, size_t count) const;

        [[nodiscard]] inline FormatedVector<Derived> format() const;
        [[nodiscard]] TransposeVector<Derived> transpose() const noexcept { return TransposeVector<Derived>(*this); }
        [[nodiscard]] ConjugateVector<Derived> conjugate() const noexcept { return ConjugateVector<Derived>(*this); }
        [[nodiscard]] HermiteVector<Derived> hermite() const noexcept { return HermiteVector<Derived>(*this); }
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }

        [[nodiscard]] inline RealType norm1() const;
        [[nodiscard]] inline RealType norm2() const;
        [[nodiscard]] inline RealType norm() const;
        [[nodiscard]] inline RealType squaredNorm() const;
        [[nodiscard]] RealType lnSquaredNorm() const;
        [[nodiscard]] inline RealType normInf() const;

        [[nodiscard]] ScalarType max() const;
        [[nodiscard]] ScalarType min() const;
        [[nodiscard]] ScalarType sum() const;
        [[nodiscard]] ScalarType prod() const;
        [[nodiscard]] bool isZeros() const;
        template<class OtherDerived>
        [[nodiscard]] inline CrossProduct<Derived, OtherDerived> crossProduct(const RValueVector<OtherDerived>& v) const noexcept;
        template<class OtherDerived>
        [[nodiscard]] ScalarType angleTo(const RValueVector<OtherDerived>& v) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> head(size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> head(size_t to) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> tail(size_t from);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> tail(size_t from) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> segment(size_t from, size_t to) const;
        [[nodiscard]] inline ReverseVector<Derived> reverse();
        [[nodiscard]] inline const ReverseVector<Derived> reverse() const;

        template<class OtherDerived>
        ReshapedVector<Derived, MatrixOption::getMajor<OtherDerived>(), OtherDerived::RowAtCompile, OtherDerived::ColumnAtCompile>
        reshape(const RValueMatrix<OtherDerived>& mat) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        ReshapedVector<Derived, MatrixOption::Column, Row, Column> reshape_col(size_t row, size_t col) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        ReshapedVector<Derived, MatrixOption::Row, Row, Column> reshape_row(size_t row, size_t col) const;
    protected:
        RValueVector() = default;
        RValueVector(const RValueVector&) = default;
        RValueVector(RValueVector&&) noexcept = default;
        /* Operators */
        RValueVector& operator=(const RValueVector&) = default;
        RValueVector& operator=(RValueVector&&) noexcept = default;
    };

    template<class VectorType1, class VectorType2>
    bool vectorNear(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2, double precision);
}

namespace Physica {
    template<class T>
    class Traits<Core::RValueVector<T>> {
    public:
        using Derived = T;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.h"
#include "RValueVectorImpl/ReverseVector.h"
#include "CrossProduct.h"
#include "InnerDot.h"
#include "VectorExpr.h"
#include "VectorConvert.h"
