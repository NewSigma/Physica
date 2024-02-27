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

namespace Physica::Core {
    template<class ScalarType, unsigned int Order>
    Differentiable<ScalarType, DiffMode::Forward, Order>::Differentiable(ScalarType value_)
            : value(std::move(value_)), grad(0) {}

    template<class ScalarType, unsigned int Order>
    Differentiable<ScalarType, DiffMode::Forward, Order>::Differentiable(ScalarType value_, ScalarType grad_)
        : value(std::move(value_)), grad(std::move(grad_)) {}

    template<class ScalarType, unsigned int Order>
    inline bool Differentiable<ScalarType, DiffMode::Forward, Order>::operator==(const This& other) const {
        return value == other.value && grad == other.grad;
    }

    template<class ScalarType, unsigned int Order>
    inline Differentiable<ScalarType, DiffMode::Forward, Order> Differentiable<ScalarType, DiffMode::Forward, Order>::operator-() const {
        return {-value, -grad};
    }

    template<class ScalarType, unsigned int Order>
    void Differentiable<ScalarType, DiffMode::Forward, Order>::swap(Differentiable& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        value.swap(obj.value);
        grad.swap(obj.grad);
    }

    template<class ScalarType, unsigned int Order>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Forward, Order>
    Differentiable<ScalarType, DiffMode::Forward, Order>::random_uniform(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_uniform(gen));
    }

    template<class ScalarType, unsigned int Order>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Forward, Order>
    Differentiable<ScalarType, DiffMode::Forward, Order>::random_normal(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_normal(gen));
    }

    template<class ScalarType, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    [[nodiscard]] inline Differentiable<ScalarType, DiffMode::Forward, Order>
    Differentiable<ScalarType, DiffMode::Forward, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Differentiable(ScalarType::random_any(dist, gen));
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType, unsigned int Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order>::Differentiable(ScalarType value)
            : Differentiable(TracerType::getInstance().pushOperation(std::move(value), ExpressionType::Set)) {}

    template<class ScalarType, unsigned int Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order>::Differentiable(
            [[maybe_unused]] ScalarType value, [[maybe_unused]] ScalarType grad) {
        throw std::runtime_error("[Error]: This function is provided for template meta programming, you should not arrive here");
    }

    template<class ScalarType, unsigned int Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order>::Differentiable(ScalarType* pValue_, ScalarType* pGrad_)
            : pValue(pValue_), pGrad(pGrad_) {}

    template<class ScalarType, unsigned int Order>
    inline bool Differentiable<ScalarType, DiffMode::Reverse, Order>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (pGrad == other.pGrad) && "[Error]: Bad scalar");
        return result;
    }

    template<class ScalarType, unsigned int Order>
    inline Differentiable<ScalarType, DiffMode::Reverse, Order> Differentiable<ScalarType, DiffMode::Reverse, Order>::operator-() const {
        auto& tracer = TracerType::getInstance();
        const auto result = tracer.pushOperation(-getValue(), ExpressionType::Minus);
        tracer.pushOperand(*this);
        return result;
    }

    template<class ScalarType, unsigned int Order>
    inline Differentiable<ScalarType, DiffMode::Reverse, Order>& Differentiable<ScalarType, DiffMode::Reverse, Order>::toAbs() {
        (*this) = abs(*this);
        return *this;
    }

    template<class ScalarType, unsigned int Order>
    inline void Differentiable<ScalarType, DiffMode::Reverse, Order>::reverse() {
        *pGrad = ScalarType(1);
        TracerType::getInstance().reverse_from(*this);
    }

    template<class ScalarType, unsigned int Order>
    inline void Differentiable<ScalarType, DiffMode::Reverse, Order>::reverse_to(Differentiable to) {
        *pGrad = ScalarType(1);
        TracerType::getInstance().reverse(*this, to);
    }

    template<class ScalarType, unsigned int Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order> Differentiable<ScalarType, DiffMode::Reverse, Order>::copy() const {
        auto& tracer = TracerType::getInstance();
        const auto result = tracer.pushOperation(getValue(), ExpressionType::Assign);
        tracer.pushOperand(*this);
        return result;
    }

    template<class ScalarType, unsigned int Order>
    inline void Differentiable<ScalarType, DiffMode::Reverse, Order>::swap(Differentiable<ScalarType, DiffMode::Reverse, Order>& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(pGrad, obj.pGrad);
    }

    template<class ScalarType, unsigned int Order>
    inline ExpressionType Differentiable<ScalarType, DiffMode::Reverse, Order>::getSource() const noexcept {
        auto& tracer = TracerType::getInstance();
        return tracer.getSource(*this);
    }

    template<class ScalarType, unsigned int Order>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Reverse, Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order>::random_uniform(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_uniform(gen));
    }

    template<class ScalarType, unsigned int Order>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Reverse, Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order>::random_normal(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_normal(gen));
    }

    template<class ScalarType, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Reverse, Order>
    Differentiable<ScalarType, DiffMode::Reverse, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Differentiable(ScalarType::random_any(dist, gen));
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type
    operator+(const Differentiable<ScalarType, Mode, 1>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode, 1>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() + s2.getValue(), s1.getGrad() + s2.getGrad());
            else
                return ResultType(s1.getValue() + s2.getValue(), s1.getGrad());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const PlainScalar value = s1.getValue() + s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExpressionType::Add);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExpressionType::Add);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, 1>& s2) {
        return s2 + s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type
    operator-(const Differentiable<ScalarType, Mode, 1>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode, 1>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() - s2.getValue(), s1.getGrad() - s2.getGrad());
            else
                return ResultType(s1.getValue() - s2.getValue(), s1.getGrad());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const PlainScalar value = s1.getValue() - s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExpressionType::Sub);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExpressionType::Sub);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, 1>& s2) {
        return -(s2 - s1.getDerived());
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type
    operator*(const Differentiable<ScalarType, Mode, 1>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode, 1>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() * s2.getValue(), s1.getGrad() * s2.getValue() + s1.getValue() * s2.getGrad());
            else
                return ResultType(s1.getValue() * s2.getValue(), s1.getGrad() * s2.getValue());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const PlainScalar value = s1.getValue() * s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExpressionType::Mul);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExpressionType::Mul);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, 1>& s2) {
        return s2 * s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type
    operator/(const Differentiable<ScalarType, Mode, 1>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode, 1>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            const auto rep = reciprocal(s2.getValue());
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() * rep, (s1.getGrad() * s2.getValue() - s1.getValue() * s2.getGrad()) * square(rep));
            else
                return ResultType(s1.getValue() * rep, s1.getGrad() * rep);
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const PlainScalar value = s1.getValue() / s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExpressionType::Div);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExpressionType::Div);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode, 1>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode, 1>, OtherScalar>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            const auto rep = reciprocal(s2.getValue());
            return ResultType(s1.getValue() * rep, -s1.getValue() * s2.getGrad() * square(rep));
        }
        else {
            static_assert(std::is_same<OtherScalar, ScalarType>::value, "[Error]: Reverse mode between different type is not supported");
            const ScalarType value = s1.getValue() / s2.getValue();
            const ResultType copy = s1.getDerived();
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Div);
            tracer.pushOperand(copy);
            tracer.pushOperand(s2.getDerived());
            return result;
        }
    }
}
