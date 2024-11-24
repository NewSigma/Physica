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

#include "Physica/CRTPBase.h"
#include "Physica/Core/MultiPrecision/Real.h"
#include "Physica/Core/MultiPrecision/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "RValueVectorImpl/RVectorBlock.h"

namespace Physica::Core {
    template<class Derived> class LValueVector;
    template<class Derived> class ContinuousVector;
    template<class Derived> class RValueMatrix;
    template<class Derived> class ContinuousMatrix;
    template<Scalar, int Option, size_t Row, size_t Col, class Allocator> class DenseMatrix;
    template<class VectorType> class TransposeVector;
    template<class VectorType> class ConjugateVector;
    template<class VectorType> class HermiteVector;
    template<class VectorType> class FormatedVector;
    template<class VectorType> class ReverseVector;
    template<class VectorType, int MatrixMajor, size_t Row, size_t Col> class ReshapedVector;
    template<Vector V1, Vector V2> class CrossProduct;

    template<class T>
    struct is_continuous {
        constexpr static bool value = std::is_base_of<ContinuousVector<T>, T>::value || std::is_base_of<ContinuousMatrix<T>, T>::value;
    };

    namespace Internal {
        template<Vector T1, Vector T2 = T1>
        class EnableSIMD {
            using ScalarType = T1::ScalarType;
            using ValueType = ScalarType::ValueType;
            constexpr static bool isSameScalar = std::is_same<ValueType, typename T2::ValueType>::value;
        public:
            constexpr static bool value = isSameScalar && BestPacket<ScalarType, T1::SizeAtCompile>::Size > 1;
        };
    }

    template<Vector T1, Vector T2 = T1>
    struct EnableMKL {
        constexpr static bool value = HasMKL() && is_continuous<T1>::value && is_continuous<T2>::value;
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
        using ScalarType = Traits<Derived>::ScalarType;
        using ValueType = ScalarType::ValueType;
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        using PacketType = BestPacket<ScalarType, SizeAtCompile>::Type;
        using ColMatrix = DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Vector, SizeAtCompile, 1, HostAllocator<ScalarType>>;
        using RowMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector, 1, SizeAtCompile, HostAllocator<ScalarType>>;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
    private:
        using RealType = ScalarType::RealType;
    public:
        ~RValueVector() = default;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& v) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        template<class AnyPacket>
        [[nodiscard]] inline AnyPacket packet(size_t index) const;
        template<class AnyPacket>
        [[nodiscard]] inline AnyPacket packetPartial(size_t index, size_t count) const;

        [[nodiscard]] inline auto format() const;
        [[nodiscard]] auto transpose() const noexcept { return TransposeVector<Derived>(Base::getDerived()); }
        [[nodiscard]] auto conjugate() const noexcept { return ConjugateVector<Derived>(Base::getDerived()); }
        [[nodiscard]] auto hermite() const noexcept { return HermiteVector<Derived>(Base::getDerived()); }
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
        [[nodiscard]] ScalarType lnSumExp() const;
        [[nodiscard]] ScalarType prod() const;
        [[nodiscard]] bool isZeros() const;
        template<Vector V>
        [[nodiscard]] inline auto crossProduct(const V& v) const noexcept;
        template<Vector V>
        [[nodiscard]] ScalarType angleTo(const V& v) const noexcept;
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

        template<Matrix M>
        ReshapedVector<Derived, MatrixOption::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>
        reshape(const M& mat) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        ReshapedVector<Derived, MatrixOption::Col, Row, Col> reshape_col(size_t row, size_t col) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        ReshapedVector<Derived, MatrixOption::Row, Row, Col> reshape_row(size_t row, size_t col) const;
    protected:
        RValueVector() = default;
        RValueVector(const RValueVector&) = default;
        RValueVector(RValueVector&&) noexcept = default;
        /* Operators */
        RValueVector& operator=(const RValueVector&) = default;
        RValueVector& operator=(RValueVector&&) noexcept = default;
    };

    template<Vector T1, Vector T2>
    bool vectorNear(const T1& v1, const T2& v2, double precision);
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
