/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Scalar/Real.h" // IWYU pragma: export
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"
#include "RValueVectorImpl/RVectorBlock.h"

namespace Physica::Core {
    template<class Derived> class LValueVector;
    template<class Derived> class ContinuousVector;
    template<class Derived> class RValueMatrix;
    template<class Derived> class ContinuousMatrix;
    template<class VectorType> class TransposeVector;
    template<class VectorType> class ConjugateVector;
    template<class VectorType> class HermiteVector;
    template<Vector V1, Vector V2> class CrossProduct;
    template<class T> class RealVector;
    template<class T> class ImagVector;
    template<class T> class SquaredNormVector;
    template<class T> class NormVector;
    template<class T> class ValueVector;
    template<class T, int GradOrder> class GradVector;

    template<class T>
    struct is_continuous {
        constexpr static bool value = std::is_base_of<ContinuousVector<T>, T>::value || std::is_base_of<ContinuousMatrix<T>, T>::value;
    };

    namespace Internal {
        template<Vector T1, Vector T2 = T1>
        class EnableSIMD {
            constexpr static size_t Size1 = T1::SizeAtCompile;
            constexpr static size_t Size2 = T2::SizeAtCompile;
        public:
            constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
            using ResultType = BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
            using PacketType = BestPacket<ResultType, SizeAtCompile>::Type;
        private:
            constexpr static bool isSameScalar = std::is_same<typename T1::ValueType, typename T2::ValueType>::value;
            constexpr static bool isBadPacket = PacketType::size() == 1;
        public:
            constexpr static bool value = isSameScalar && !isBadPacket;
        };

        template<Vector T1, Vector T2 = T1>
        struct EnableMKL {
            constexpr static bool value = HasMKL()
                                       && is_continuous<T1>::value
                                       && is_continuous<T2>::value
                                       && !Diffable<T1>
                                       && !Diffable<T2>
                                       && (EnableSIMD<T1, T2>::SizeAtCompile == Dynamic);
        };
    }
    /**
     * \class RValueVector: The base class for all vectors.
     */
    template<class Derived>
    class RValueVector : public CRTPBase<RValueVector<Derived>> {
        using This = RValueVector<Derived>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
        using ValueType = ScalarType::ValueType;
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        using PacketType = BestPacket<ScalarType, SizeAtCompile>::Type;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
    private:
        template<size_t Length>
        using BlockType = RVectorBlock<Derived, Length>;
        using RealType = ScalarType::RealType;
    public:
        ~RValueVector() = default;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& v) const;

        [[nodiscard]] auto calc(size_t index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] auto calc_value(size_t index) const { return Base::getDerived().calc_value(index); }
        template<class AnyPacket>
        [[nodiscard]] inline AnyPacket packet(size_t index) const;
        template<class AnyPacket>
        [[nodiscard]] inline AnyPacket packetPartial(size_t index, size_t count) const;

        [[nodiscard]] inline auto format() const;
        [[nodiscard]] auto transpose() const noexcept;
        [[nodiscard]] auto conjugate() const noexcept;
        [[nodiscard]] auto hermite() const noexcept;

        [[nodiscard]] inline CoDiff<RealType> norm1() const;
        [[nodiscard]] inline CoDiff<RealType> norm2() const;
        [[nodiscard]] inline CoDiff<RealType> norm() const;
        [[nodiscard]] inline CoDiff<RealType> squaredNorm() const;
        [[nodiscard]] inline RealType lnSquaredNorm() const;
        [[nodiscard]] inline RealType normInf() const;

        [[nodiscard]] CoDiff<ScalarType> max() const;
        [[nodiscard]] CoDiff<ScalarType> min() const;
        [[nodiscard]] CoDiff<ScalarType> sum() const;
        [[nodiscard]] CoDiff<ScalarType> lnSumExp() const;
        [[nodiscard]] CoDiff<ScalarType> crossEntropy(size_t index) const;
        [[nodiscard]] CoDiff<ScalarType> lnSoftmax(size_t index) const;
        [[nodiscard]] CoDiff<ScalarType> softmax(size_t index) const;
        [[nodiscard]] CoDiff<ScalarType> prod() const;
        [[nodiscard]] bool isZeros() const;
        template<Vector V>
        [[nodiscard]] inline auto crossProduct(const V& v) const noexcept;
        template<Vector V>
        [[nodiscard]] ScalarType angleTo(const V& v) const noexcept;

        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto segment(size_t from, size_t to) const noexcept;
        [[nodiscard]] inline auto reversal() noexcept;
        [[nodiscard]] inline const auto reversal() const noexcept;

        template<Matrix M>
        auto reshape(const M& mat) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        auto reshape_col(size_t row, size_t col) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        auto reshape_row(size_t row, size_t col) const noexcept;

        auto reals() const noexcept;
        auto imags() const noexcept;
        auto squaredNorms() const noexcept;
        auto norms() const noexcept;
        auto values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
    protected:
        RValueVector() = default;
        RValueVector(const RValueVector&) = default;
        RValueVector(RValueVector&&) noexcept = default;
        /* Operators */
        RValueVector& operator=(const RValueVector&) = default;
        RValueVector& operator=(RValueVector&&) noexcept = default;
        /* Operations */
        template<int GradOrder>
        auto grads_impl() const noexcept;
    };

    template<Vector T1, Vector T2>
    bool vectorNear(const T1& v1, const T2& v2, double precision);
}

namespace Physica {
    template<class T>
    class Traits<RValueVector<T>> {
    public:
        using Derived = T;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.h"
#include "RValueVectorImpl/ReversalVector.h"
#include "RValueVectorImpl/CrossProduct.h"
#include "RValueVectorImpl/VectorConvert.h"
#include "InnerDot.h"
#include "VectorExpr.h"
