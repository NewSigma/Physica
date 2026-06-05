/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Parallel/Parallel.h"
#include "RValueVectorImpl/RVectorBlock.h"

namespace Physica {
    template<class Derived> class LValueVector;
    template<class Derived> class CompactVector;
    template<class Derived> class RValueMatrix;
    template<class Derived> class CompactMatrix;
    template<Scalar, size_t Length, class Allocator> class DenseVector;
    template<Vector, int MatrixMajor, size_t Row, size_t Col> class RValueReshapedVector;
    template<class> class Transpose;
    template<class> class Conjugate;
    template<class> class Hermite;
    template<Vector, Vector> class CrossProduct;
    template<class T> class RealVectorR;
    template<class T> class ImagVector;
    template<class T> class SquaredNormVector;
    template<class T> class NormVector;
    template<class T> class ValueVector;
    template<class T, int GradOrder> class GradVector;
    template<class T, int MaskOrder> class GradMaskVector;
    template<Vector, Matrix> class GEVM;
    template<Matrix, Vector> class GEMV;

    namespace Internal {
        /**
         * LAPACK has requirements on its parameters; we can check them using \class EnableLAPACK
         */
        template<class T1, class T2 = T1> class EnableLAPACK;

        template<Vector V1, Vector V2 = V1>
        class EnableSIMD {
            using U1 = std::remove_cvref<V1>::type;
            using U2 = std::remove_cvref<V2>::type;
            using T1 = typename U1::ScalarType;
            using T2 = typename U2::ScalarType;
        public:
            constexpr static size_t SizeAtCompile = std::max(U1::getSizeAtCompile(), U2::getSizeAtCompile());
            using ResultType = BinaryScalarOpRtnTy<T1, T2>::Type;
            using PacketType = BestPacket<ResultType, SizeAtCompile>::Type;

            constexpr static bool value = []() consteval static noexcept {
                constexpr bool isSameScalar = std::same_as<typename T1::ValueType, typename T2::ValueType>;
                constexpr bool isCUDA = DeviceObj<V1> || DeviceObj<V2>;
                constexpr bool isFloat16 = ResultType::Prec == Float16;
                // Only use FP16 SIMD for device:
                // 1. Other packet types do not work for CUDA
                // 2. Old processors do not have FP16 support
                constexpr bool UsePacketFP16 = isFloat16 == isCUDA;
                constexpr bool FastPacket = U1::isFastPacket() && U2::isFastPacket();
                return isSameScalar && !Scalar<PacketType> && UsePacketFP16 && FastPacket;
            }();
        };

        template<Vector V1, Vector V2>
        class EnableLAPACK<V1, V2> {
            using U1 = std::remove_cvref<V1>::type;
            using U2 = std::remove_cvref<V2>::type;
            using T = U1::ScalarType;
        public:
            constexpr static bool value = std::same_as<T, typename U2::ScalarType>
                                       && (T::Prec == Float16 || T::Prec == Float32 || T::Prec == Float64)
                                       && U1::isCompact()
                                       && U2::isCompact()
                                       && !Diffable<U1>
                                       && (EnableSIMD<U1, U2>::SizeAtCompile == Dynamic);
        };
    }
    /**
     * \class RValueVector: The base class for all vectors.
     */
    template<class Derived>
    class RValueVector : public CRTPBase<RValueVector<Derived>> {
        static_assert(!DeviceObj<Derived>, "[Error]: device_obj<> must be outside RValueVector<>");
        using This = RValueVector<Derived>;
        using Base = CRTPBase<This>;

        template<Vector> class View;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
        using Tcv = Tc::ValueType;
    public:
        ~RValueVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] bool operator!=(const Vector auto& other) const noexcept;
        [[nodiscard, gnu::always_inline]] auto operator*(this auto&&, Scalar auto&& x) noexcept;
        [[nodiscard]] auto operator*(this auto&&, Matrix auto&& m) noexcept;
        [[nodiscard, gnu::always_inline]] auto operator-(this auto&&) noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Vector auto&& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto&& v) const noexcept;
        void assert_assign(const Vector auto& source) const noexcept;
        void assert_assign_lapack(const Vector auto& source) const noexcept;

        [[nodiscard]] decltype(auto) calc(size_t index) const noexcept;
        [[nodiscard]] decltype(auto) calc_value(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
        [[nodiscard]] constexpr auto view(this auto&&) noexcept;
        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;
        void resize(const Vector auto& x);
        void resize(size_t length);

        template<size_t Length = Dynamic>
        [[nodiscard]] auto head(this auto&&, size_t to = Length) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto segment(this auto&&, size_t from, size_t to) noexcept;
        template<int Major, size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape(size_t row, size_t col) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_row(size_t row, size_t col) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_col(size_t row, size_t col) const noexcept;
        [[nodiscard]] auto reshape_like(const Matrix auto& mat) const noexcept;
        [[nodiscard]] auto reversal() noexcept;
        [[nodiscard]] const auto reversal() const noexcept;
        [[nodiscard]] auto transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) conjugate(this auto&&) noexcept;
        [[nodiscard]] auto hermite(this auto&&) noexcept;
        [[nodiscard]] auto format() const;

        [[nodiscard]] CoDiff<Tr> norm1() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm2() const noexcept;
        [[nodiscard]] auto norm() const noexcept;
        [[nodiscard]] CoDiff<Tr> squaredNorm() const noexcept;
        [[nodiscard]] Tr lnSquaredNorm() const;
        [[nodiscard]] Tr normInf() const;
        [[nodiscard]] size_t argmax() const noexcept;
        [[nodiscard]] size_t argmin() const noexcept;
        [[nodiscard]] CoDiff<T> max() const noexcept;
        [[nodiscard]] CoDiff<T> min() const noexcept;
        [[nodiscard]] CoDiff<T> sum(this auto&&) noexcept;
        [[nodiscard]] CoDiff<T> mean() const noexcept;
        [[nodiscard]] CoDiff<T> mean_stable() const noexcept;
        [[nodiscard]] CoDiff<T> variance() const noexcept;
        [[nodiscard]] T variance(const T& prior_mean) const;
        [[nodiscard]] T variance_stable() const;
        [[nodiscard]] T deviation() const;
        [[nodiscard]] T deviation(const T& prior_mean) const;
        [[nodiscard]] T deviation_stable() const;
        [[nodiscard]] CoDiff<T> lnSumExp(this const auto&) noexcept;
        [[nodiscard]] CoDiff<T> crossEntropy(size_t index) const noexcept;
        [[nodiscard]] CoDiff<T> lnSoftmax(size_t index) const noexcept;
        [[nodiscard]] CoDiff<T> softmax(size_t index) const noexcept;
        [[nodiscard]] CoDiff<T> prod() const noexcept;
        [[nodiscard]] auto cross(const Vector auto& v) const noexcept;
        Tr householder(Vector auto& __restrict target) const __restrict;

        [[nodiscard]] decltype(auto) reals(this auto&&) noexcept;
        [[nodiscard]] auto imags(this auto&&) noexcept;
        [[nodiscard]] auto squaredNorms(this auto&&) noexcept;
        [[nodiscard]] auto norms(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        template<int MaskOrder>
        [[nodiscard]] decltype(auto) grads_mask(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        [[nodiscard]] bool empty() const noexcept { return getLength() == 0; }
        [[nodiscard]] bool isZero() const;
        [[nodiscard]] bool isFinite() const;
        [[nodiscard]] bool isSubNormal() const;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isForwardDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isReverseDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isLValueVector() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isSparse() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastPacket() noexcept;
        [[nodiscard]] __host__ __device__ consteval static auto getSizeAtCompile() noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Scalar auto& source) noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Vector auto& source) noexcept;
    protected:
        RValueVector() = default;
        RValueVector(const This&) = default;
        RValueVector(This&&) noexcept = default;
        /* Static members */
        [[nodiscard]] consteval static int calcBlockingSize(int CacheSize) noexcept;
    private:
        template<ExecutePolicy P>
        void assign_for(Vector auto& v) const noexcept;
        template<Vector V, ExecutePolicy P, size_t Size>
        void assign_simd(V& v) const noexcept;

        template<ExecutePolicy P>
        void assign_add_for(Vector auto& v) const noexcept;
        template<Vector V, size_t Length>
        void assign_add_simd(V& v) const noexcept;
    };

    template<Vector V>
    auto covariance(const V& x, const V& y) -> V::ScalarType {
        assert(x.getLength() == y.getLength());
        using T = V::ScalarType;
        const T x_mean = x.mean();
        const T y_mean = y.mean();
        return hadamard((x - x_mean), (y - y_mean)).sum() / T(x.getLength() - 1);
    }

    template<Vector V1, Vector V2>
    bool vectorNear(const V1& v1, const V2& v2, double precision) noexcept {
        using T = Internal::BinaryScalarOpRtnTy<typename V1::ScalarType, typename V2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        for (size_t i = 0; i < v1.getLength(); ++i)
            if (!scalarNear(T(v1.calc(i)), T(v2.calc(i)), precision))
                return false;
        return true;
    }

    template<Vector V1, Vector V2>
    bool vectorNear(const V1& v1, const V2& v2, uint64_t ulp) noexcept {
        using T = Internal::BinaryScalarOpRtnTy<typename V1::ScalarType, typename V2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        for (size_t i = 0; i < v1.getLength(); ++i)
            if (!scalarNear(T(v1.calc(i)), T(v2.calc(i)), ulp))
                return false;
        return true;
    }

    std::ostream& operator<<(std::ostream& os, const Vector auto& v) {
        return os << std::format("{}", v.format());
    }

    [[nodiscard, gnu::always_inline]] auto operator*(Scalar auto&& x, Vector auto&& v) noexcept {
        return std::forward<decltype(v)>(v) * std::forward<decltype(x)>(x);
    }
}

namespace Physica {
    template<class T>
    class Traits<RValueVector<T>> {
    public:
        using Derived = T;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.h"
#include "RValueVectorImpl/View.h"
#include "RValueVectorImpl/ReversalVector.h"
#include "RValueVectorImpl/Conjugate.h"
#include "RValueVectorImpl/CrossProduct.h"
#include "RValueVectorImpl/Dot.h"
#include "RValueVectorImpl/VectorConvert/RealVector.h"
#include "RValueVectorImpl/VectorConvert/SquaredNormVector.h"
#include "RValueVectorImpl/VectorConvert/VectorConvert.h"
#include "VectorExpr.h"
