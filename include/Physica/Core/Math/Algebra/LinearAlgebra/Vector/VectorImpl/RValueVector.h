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
#include "Physica/Core/Scalar/Complex.h" // IWYU pragma: export
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"
#include "RValueVectorImpl/RVectorBlock.h"

namespace Physica {
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
    class is_continuous {
        using U = std::remove_cvref<T>::type;
    public:
        constexpr static bool value = std::is_base_of<ContinuousVector<U>, U>::value || std::is_base_of<ContinuousMatrix<U>, U>::value;
    };

    namespace Internal {
        template<Vector T1, Vector T2 = T1>
        class EnableSIMD {
            constexpr static size_t Size1 = T1::SizeAtCompile;
            constexpr static size_t Size2 = T2::SizeAtCompile;
            using U1 = typename T1::ScalarType;
            using U2 = typename T2::ScalarType;
        public:
            constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
            using ResultType = BinaryScalarOpRtnTy<U1, U2>::Type;
            using PacketType = BestPacket<ResultType, SizeAtCompile>::Type;
        private:
            constexpr static bool isSameScalar = std::same_as<typename U1::ValueType, typename U2::ValueType>;
            constexpr static bool isBadPacket = PacketType::size() == 1;
            constexpr static bool isCUDA = CUDA<T1> || CUDA<T2>;
            constexpr static bool isFloat16 = ResultType::Option == Float16;
        public:
            constexpr static bool value = (isCUDA == isFloat16) && isSameScalar && !isBadPacket;
        };

        template<Vector T1, Vector T2 = T1>
        class EnableMKL {
            using U1 = std::remove_cvref<T1>::type;
            using U2 = std::remove_cvref<T2>::type;
            using ScalarType1 = U1::ScalarType;
            using ScalarType2 = U2::ScalarType;
        public:
            constexpr static bool value = HasMKL()
                                       && std::same_as<ScalarType1, ScalarType2>
                                       && (ScalarType1::Option == Float32 || ScalarType1::Option == Float64)
                                       && is_continuous<U1>::value
                                       && is_continuous<U2>::value
                                       && !Diffable<U1>
                                       && !Diffable<U2>
                                       && (EnableSIMD<U1, U2>::SizeAtCompile == Dynamic);
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
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        using PacketType = BestPacket<ScalarType, SizeAtCompile>::Type;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isDiffable = ScalarType::isDiffable;
        constexpr static bool isComplex = ScalarType::isComplex;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
    private:
        template<size_t Length>
        using BlockType = RVectorBlock<Derived, Length>;
        using RealsRtnTy = std::conditional<isComplex, RealVector<Derived>, Derived&>::type;
        using ValuesRtnTy = std::conditional<isDiffable, ValueVector<Derived>, Derived&>::type;
    public:
        ~RValueVector() = default;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& v) const;
        template<Vector V, class Executor = SeqExecutor>
        inline void assign_add(V& v) const;

        [[nodiscard]] auto calc(size_t index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] auto calc_value(size_t index) const { return Base::getDerived().calc_value(index); }
        template<Packet Pack>
        [[nodiscard]] inline Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] inline Pack packetPartial(size_t index, size_t count) const;
        template<Vector V1, Vector V2>
        void reverse(const V1& y, const V2& grad) const noexcept requires(isReverseDiff);

        [[nodiscard]] inline auto format() const;
        [[nodiscard]] auto transpose() const noexcept;
        [[nodiscard]] auto conjugate() const noexcept;
        [[nodiscard]] auto hermite() const noexcept;

        [[nodiscard]] inline CoDiff<Tr> norm1() const;
        [[nodiscard]] inline CoDiff<Tr> norm2() const;
        [[nodiscard]] inline CoDiff<Tr> norm() const;
        [[nodiscard]] inline CoDiff<Tr> squaredNorm() const;
        [[nodiscard]] inline Tr lnSquaredNorm() const;
        [[nodiscard]] inline Tr normInf() const;

        [[nodiscard]] CoDiff<T> max() const;
        [[nodiscard]] CoDiff<T> min() const;
        [[nodiscard]] CoDiff<T> sum() const;
        [[nodiscard]] CoDiff<T> mean() const;
        [[nodiscard]] CoDiff<T> variance() const;
        [[nodiscard]] T variance(const T& prior_mean) const;
        [[nodiscard]] T deviation() const;
        [[nodiscard]] T deviation(const T& prior_mean) const;
        [[nodiscard]] CoDiff<T> lnSumExp() const;
        [[nodiscard]] CoDiff<T> crossEntropy(size_t index) const;
        [[nodiscard]] CoDiff<T> lnSoftmax(size_t index) const;
        [[nodiscard]] CoDiff<T> softmax(size_t index) const;
        [[nodiscard]] CoDiff<T> prod() const;
        [[nodiscard]] bool isZeros() const;
        [[nodiscard]] bool isFinite() const;
        template<Vector V>
        [[nodiscard]] inline auto crossProduct(const V& v) const noexcept;
        template<Vector V>
        [[nodiscard]] T angleTo(const V& v) const noexcept;

        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto head(size_t to) & noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto head(size_t to) const& noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto tail(size_t from) & noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto tail(size_t from) const& noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto segment(size_t from, size_t to) & noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto segment(size_t from, size_t to) const& noexcept;
        [[nodiscard]] inline auto reversal() noexcept;
        [[nodiscard]] inline const auto reversal() const noexcept;

        template<Matrix M>
        auto reshape(const M& mat) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        auto reshape_col(size_t row, size_t col) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        auto reshape_row(size_t row, size_t col) const noexcept;

        [[nodiscard]] RealsRtnTy reals() const noexcept;
        [[nodiscard]] auto imags() const noexcept;
        [[nodiscard]] auto squaredNorms() const noexcept;
        [[nodiscard]] auto norms() const noexcept;
        [[nodiscard]] ValuesRtnTy values() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        /* Static members */
        template<Vector V>
        __host__ __device__ static void assign_check(const V& target) noexcept;
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
    private:
        template<Vector V, class Executor>
        inline void assign_for(V& v) const;
        template<Vector V, class Executor, size_t Size>
        inline void assign_simd(V& v) const;

        template<Vector V, class Executor>
        inline void assign_add_for(V& v) const;
        template<Vector V, size_t Size>
        inline void assign_add_simd(V& v) const;
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
#include "RValueVectorImpl/InnerDot.h"
#include "VectorExpr.h"
