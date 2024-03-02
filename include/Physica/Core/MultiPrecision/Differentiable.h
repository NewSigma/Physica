/*
 * Copyright 2023-2024 WeiBo He.
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
#include "DifferentiableImpl/DiffTracer.h"

namespace Physica::Core {
    namespace Internal {
        template<class T, DiffMode M, unsigned int Order_>
        class Traits<Differentiable<T, M, Order_>> {
            static_assert(!T::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
            static_assert(Order_ > 0, "[Error]: Use plain type instead of 0 order differentiable");
            using RealT = typename T::RealType;
            using ComplexT = typename T::ComplexType;
        public:
            using PlainScalar = T;
            constexpr static DiffMode Mode = M;
            constexpr static unsigned int Order = Order_;
            using ScalarType = Differentiable<T, M, Order>;
            using RealType = Differentiable<RealT, M, Order>;
            using ComplexType = Differentiable<ComplexT, M, Order>;
            using TrivialType = typename T::TrivialType;
            constexpr static ScalarOption Option = T::Option;
            constexpr static bool isComplex = T::isComplex;
            constexpr static bool isDifferentiable = true;
            constexpr static bool isForwardDiff = Mode == DiffMode::Forward;
            constexpr static bool isReverseDiff = Mode == DiffMode::Reverse;
            /* SIMD */
            using BoolSIMDType = BoolSIMD<ScalarType, 1>;
        };

        template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
        class BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar> {
        public:
            using Type = Differentiable<typename BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type, Mode, Order>;
        };

        template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
        class BinaryScalarOpReturnType<OtherScalar, Differentiable<ScalarType, Mode, Order>> {
        public:
            using Type = typename BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type;
        };

        template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar, DiffMode OtherMode, unsigned int OtherOrder>
        class BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, Differentiable<OtherScalar, OtherMode, OtherOrder>> {
            static_assert(Mode == OtherMode, "[Error]: Operation between differentiable scalars with different mode is not well defined");
            static_assert(Order == OtherOrder, "[Error]: Operation between differentiable scalars with different order is not well defined");
        public:
            using Type = Differentiable<typename BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type, Mode, Order>;
        };
    }
    /**
     * \class Differentiable provides auto differential support for scalars
     */
    template<class ScalarType, unsigned int Order>
    class Differentiable<ScalarType, DiffMode::Forward, Order> : public ScalarBase<Differentiable<ScalarType, DiffMode::Forward, Order>> {
        static_assert(Order == 1, "[Error]: High order autodiff is not implemented");
        using This = Differentiable<ScalarType, DiffMode::Forward, Order>;

        ScalarType value;
        ScalarType grad;
    public:
        Differentiable() = default;
        Differentiable(double d) : This(ScalarType(d)) {}
        Differentiable(ScalarType value_);
        Differentiable(ScalarType value_, ScalarType grad_);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(Differentiable obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] explicit operator float() const { return float(value); }
        [[nodiscard]] explicit operator double() const { return double(value); }
        [[nodiscard]] inline bool operator==(const This& other) const;
        [[nodiscard]] inline Differentiable operator-() const;
        /* Operations */
        void swap(Differentiable& __restrict obj) noexcept;
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
        [[nodiscard]] inline static Differentiable random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_any(Distribution& dist, RandomGenerator& gen);
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType, unsigned int Order>
    class Differentiable<ScalarType, DiffMode::Reverse, Order> : public ScalarBase<Differentiable<ScalarType, DiffMode::Reverse, Order>> {
        using This = Differentiable<ScalarType, DiffMode::Reverse, Order>;
        using ReducedType = Differentiable<ScalarType, DiffMode::Reverse, Order - 1>;
        using ValueType = ScalarType* __restrict;
        using GradType = typename std::conditional<Order == 1, ValueType, ReducedType>::type;
        template<unsigned int GradOrder>
        using GradReturnType = typename std::conditional<Order == GradOrder
                                                        , ScalarType&
                                                        , Differentiable<ScalarType, DiffMode::Reverse, Order - GradOrder>&>::type;
    public:
        using device_obj_type = device_obj<This>;
        using TracerType = DiffTracer<ScalarType, Order>;
    private:
        ValueType pValue;
        GradType grad;
    public:
        Differentiable() = default;
        Differentiable(double d) : This(ScalarType(d)) {}
        Differentiable(ScalarType value);
        Differentiable(ScalarType value, ScalarType grad);
        Differentiable(ValueType pValue_, GradType grad_);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(const Differentiable&) = default;
        Differentiable& operator=(Differentiable&&) noexcept = default;
        [[nodiscard]] explicit operator float() const { return float(getValue()); }
        [[nodiscard]] explicit operator double() const { return double(getValue()); }
        template<unsigned int KeepDeep>
        [[nodiscard]] explicit operator Differentiable<ScalarType, DiffMode::Reverse, KeepDeep>() const;
        [[nodiscard]] inline bool operator==(const This& other) const;
        [[nodiscard]] inline Differentiable operator-() const;
        /* Operations */
        inline Differentiable& toAbs();
        inline void reverse();
        inline void reverse_to(Differentiable to);
        inline void zero_grad();
        [[nodiscard]] inline This copy() const;
        inline void swap(Differentiable& __restrict obj) noexcept;
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
        [[nodiscard]] inline ExpressionType getSource() const noexcept;
        [[nodiscard]] inline This* getFirstOperand() const;
        /* Setters */
        void setValue(const ScalarType& x) { *pValue = x; }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_any(Distribution& dist, RandomGenerator& gen);
    private:
        /* Friends */
        friend class device_obj<This>;
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type
    operator+(const Differentiable<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, Order>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type
    operator-(const Differentiable<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, Order>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type
    operator*(const Differentiable<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, Order>& s2);
    
    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type
    operator/(const Differentiable<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, Order>& s2);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    inline std::ostream& operator<<(std::ostream& os, const Differentiable<ScalarType, Mode, Order>& obj) {
        return os << obj.getValue();
    }
}

namespace std {
    template<class ScalarType, Physica::Core::DiffMode Mode, unsigned int Order>
    struct numeric_limits<Physica::Core::Differentiable<ScalarType, Mode, Order>> : public numeric_limits<ScalarType> {};
}

#include "DifferentiableImpl/DifferentiableImpl.h"
#include "DifferentiableImpl/ElementaryFunction.h"
#include "DifferentiableImpl/ProbabilityFunction.h"
