/*
 * Copyright 2023-2024 Weibo He.
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

#include "Scalar.h"
#include "DiffImpl/DiffTracer.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
        class BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar> {
        public:
            using Type = Diff<typename BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type, Mode, Order>;
        };

        template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
        class BinaryScalarOpReturnType<OtherScalar, Diff<ScalarType, Mode, Order>> {
        public:
            using Type = typename BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type;
        };

        template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar, DiffMode OtherMode, unsigned int OtherOrder>
        class BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, Diff<OtherScalar, OtherMode, OtherOrder>> {
            static_assert(Mode == OtherMode, "[Error]: Operation between differentiable scalars with different mode is not well defined");
            static_assert(Order == OtherOrder, "[Error]: Operation between differentiable scalars with different order is not well defined");
        public:
            using Type = Diff<typename BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type, Mode, Order>;
        };
    }
    /**
     * \class Diff provides auto differential support for scalars
     */
    template<class ScalarType, unsigned int Order>
    class Diff<ScalarType, DiffMode::Forward, Order> : public ScalarBase<Diff<ScalarType, DiffMode::Forward, Order>> {
        static_assert(Order == 1, "[Error]: High order autodiff is not implemented");
        using This = Diff<ScalarType, DiffMode::Forward, Order>;

        ScalarType value;
        ScalarType grad;
    public:
        Diff() = default;
        Diff(double d) : This(ScalarType(d)) {}
        Diff(ScalarType value_);
        Diff(ScalarType value_, ScalarType grad_);
        Diff(const This&) = default;
        Diff(This&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] explicit operator float() const { return float(value); }
        [[nodiscard]] explicit operator double() const { return double(value); }
        [[nodiscard]] inline bool operator==(const This& other) const;
        [[nodiscard]] inline This operator-() const;
        /* Operations */
        [[nodiscard]] This conjugate() const { noImpl(); }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType& getValue() noexcept { return value; }
        [[nodiscard]] ScalarType& getGrad() noexcept { return grad; }
        [[nodiscard]] const ScalarType& getValue() const noexcept { return value; }
        [[nodiscard]] const ScalarType& getGrad() const noexcept { return grad; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return value.isZero(); }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return value.isPositive(); }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return value.isNegative(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Diff random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static Diff random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Diff random_any(Distribution& dist, RandomGenerator& gen);
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType, unsigned int Order>
    class Diff<ScalarType, DiffMode::Reverse, Order> : public ScalarBase<Diff<ScalarType, DiffMode::Reverse, Order>> {
        using This = Diff<ScalarType, DiffMode::Reverse, Order>;
        using ReducedType = Diff<ScalarType, DiffMode::Reverse, Order - 1>;
        using ValueType = ScalarType* __restrict;
        using GradType = typename std::conditional<Order == 1, ValueType, ReducedType>::type;
        template<unsigned int GradOrder>
        using GradReturnType = typename std::conditional<Order == GradOrder
                                                        , ScalarType&
                                                        , Diff<ScalarType, DiffMode::Reverse, Order - GradOrder>&>::type;
    public:
        using device_obj_type = device_obj<This>;
        using TracerType = DiffTracer<ScalarType, Order>;
    private:
        ValueType pValue;
        GradType grad;
    public:
        Diff() = default;
        Diff(double d) : This(ScalarType(d)) {}
        Diff(ScalarType value);
        Diff(ScalarType value, ScalarType grad);
        Diff(ValueType pValue_, GradType grad_);
        Diff(const This&) = default;
        Diff(This&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        inline This& operator=(std::nullptr_t) noexcept;
        [[nodiscard]] explicit operator float() const { return float(getValue()); }
        [[nodiscard]] explicit operator double() const { return double(getValue()); }
        template<unsigned int KeepDeep>
        [[nodiscard]] explicit operator Diff<ScalarType, DiffMode::Reverse, KeepDeep>() const;
        [[nodiscard]] inline bool operator==(const This& other) const;
        [[nodiscard]] inline bool operator==(std::nullptr_t) const noexcept;
        [[nodiscard]] inline bool operator!=(std::nullptr_t) const noexcept { return !((*this) == nullptr); }
        [[nodiscard]] inline Diff operator-() const;
        /* Operations */
        [[nodiscard]] This conjugate() const { noImpl(); }
        inline void reverse();
        inline void reverse_to(Diff to);
        inline void zero_grad();
        [[nodiscard]] inline This copy() const;
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ inline ScalarType* value_ptr() const noexcept;
        [[nodiscard]] __host__ __device__ inline ScalarType* grad_ptr() const noexcept;
        [[nodiscard]] inline ScalarType& getValue() noexcept;
        template<unsigned int GradOrder = 1>
        [[nodiscard]] inline GradReturnType<GradOrder>& getGrad() noexcept;
        [[nodiscard]] inline const ScalarType& getValue() const noexcept;
        template<unsigned int GradOrder = 1>
        [[nodiscard]] inline const GradReturnType<GradOrder>& getGrad() const noexcept;
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return getValue().isZero(); }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return getValue().isPositive(); }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return getValue().isNegative(); }
        [[nodiscard]] inline ExprType getSource() const noexcept;
        [[nodiscard]] inline This* getFirstOperand() const;
        /* Setters */
        void setValue(const ScalarType& x) { *pValue = x; }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Diff random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static Diff random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Diff random_any(Distribution& dist, RandomGenerator& gen);
    private:
        /* Friends */
        friend class device_obj<This>;
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type
    operator+(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type
    operator-(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type
    operator*(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2);
    
    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type
    operator/(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    inline std::ostream& operator<<(std::ostream& os, const Diff<ScalarType, Mode, Order>& obj) {
        return os << obj.getValue();
    }
}

namespace Physica {
    using namespace Core;

    template<class T, DiffMode M, unsigned int Order_>
    class Traits<Diff<T, M, Order_>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
        static_assert(Order_ > 0, "[Error]: Use plain type instead of 0 order differentiable");
        using RealT = typename T::RealType;
        using ComplexT = typename T::ComplexType;
    public:
        using PlainScalar = T;
        constexpr static DiffMode Mode = M;
        constexpr static unsigned int Order = Order_;
        using ScalarType = Diff<T, M, Order>;
        using RealType = Diff<RealT, M, Order>;
        using ComplexType = Diff<ComplexT, M, Order>;
        using TrivialType = typename T::TrivialType;
        constexpr static ScalarOption Option = T::Option;
        constexpr static bool isComplex = T::isComplex;
        constexpr static bool isDifferentiable = true;
        constexpr static bool isForwardDiff = Mode == DiffMode::Forward;
        constexpr static bool isReverseDiff = Mode == DiffMode::Reverse;
        /* SIMD */
        using BoolSIMDType = BoolSIMD<ScalarType, 1>;
    };
}

namespace std {
    template<class ScalarType, Physica::Core::DiffMode Mode, unsigned int Order>
    struct numeric_limits<Physica::Core::Diff<ScalarType, Mode, Order>> : public numeric_limits<ScalarType> {};
}

#include "DiffImpl/DiffImpl.h"
#include "DiffImpl/Math.h"
