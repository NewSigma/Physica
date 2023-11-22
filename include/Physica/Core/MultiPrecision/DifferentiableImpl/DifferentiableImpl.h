/*
 * Copyright 2023 WeiBo He.
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
    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Forward>::Differentiable(ScalarType value_)
            : value(std::move(value_)), tangent(0) {}

    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Forward>::Differentiable(ScalarType value_, ScalarType tangent_)
        : value(std::move(value_)), tangent(std::move(tangent_)) {}

    template<class ScalarType>
    inline bool Differentiable<ScalarType, DiffMode::Forward>::operator==(const This& other) const {
        return value == other.value && tangent == other.tangent;
    }

    template<class ScalarType>
    inline Differentiable<ScalarType, DiffMode::Forward> Differentiable<ScalarType, DiffMode::Forward>::operator-() const {
        return {-value, -tangent};
    }

    template<class ScalarType>
    void Differentiable<ScalarType, DiffMode::Forward>::swap(Differentiable& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        value.swap(obj.value);
        tangent.swap(obj.tangent);
    }

    template<class ScalarType>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Forward>
    Differentiable<ScalarType, DiffMode::Forward>::random_uniform(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_uniform(gen));
    }

    template<class ScalarType>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Forward>
    Differentiable<ScalarType, DiffMode::Forward>::random_normal(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_normal(gen));
    }

    template<class ScalarType>
    template<class Distribution, class RandomGenerator>
    [[nodiscard]] inline Differentiable<ScalarType, DiffMode::Forward>
    Differentiable<ScalarType, DiffMode::Forward>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Differentiable(ScalarType::random_any(dist, gen));
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Reverse>::Differentiable(ScalarType value)
            : Differentiable(DiffTracerType::getInstance().pushOperation(std::move(value), ExpressionType::Set)) {}

    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Reverse>::Differentiable(
            [[maybe_unused]] ScalarType value, [[maybe_unused]] ScalarType tangent) {
        throw std::runtime_error("[Error]: This function is provided for template meta programming, you should not arrive here");
    }

    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Reverse>::Differentiable(ScalarType* pValue_, ScalarType* pTangent_)
            : pValue(pValue_), pTangent(pTangent_) {}

    template<class ScalarType>
    inline bool Differentiable<ScalarType, DiffMode::Reverse>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (pTangent == other.pTangent) && "[Error]: Bad scalar");
        return result;
    }

    template<class ScalarType>
    inline Differentiable<ScalarType, DiffMode::Reverse> Differentiable<ScalarType, DiffMode::Reverse>::operator-() const {
        auto& tracer = DiffTracerType::getInstance();
        const auto result = tracer.pushOperation(-getValue(), ExpressionType::Minus);
        tracer.pushOperand(*this);
        return result;
    }

    template<class ScalarType>
    inline void Differentiable<ScalarType, DiffMode::Reverse>::reverse() {
        auto& tracer = DiffTracer<ScalarType>::getInstance();
        *pTangent = ScalarType(1);
        tracer.reverse(*this);
    }

    template<class ScalarType>
    inline void Differentiable<ScalarType, DiffMode::Reverse>::reverse(Differentiable to) {
        auto& tracer = DiffTracer<ScalarType>::getInstance();
        *pTangent = ScalarType(1);
        tracer.reverse(*this, to);
    }

    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Reverse> Differentiable<ScalarType, DiffMode::Reverse>::copy() const {
        auto& tracer = DiffTracer<ScalarType>::getInstance();
        const auto result = tracer.pushOperation(getValue(), ExpressionType::Assign);
        tracer.pushOperand(*this);
        return result;
    }

    template<class ScalarType>
    inline void Differentiable<ScalarType, DiffMode::Reverse>::swap(Differentiable<ScalarType, DiffMode::Reverse>& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(pTangent, obj.pTangent);
    }

    template<class ScalarType>
    inline ExpressionType Differentiable<ScalarType, DiffMode::Reverse>::getSource() const noexcept {
        auto& tracer = DiffTracer<ScalarType>::getInstance();
        return tracer.getSource(*this);
    }

    template<class ScalarType>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Reverse>
    Differentiable<ScalarType, DiffMode::Reverse>::random_uniform(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_uniform(gen));
    }

    template<class ScalarType>
    template<class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Reverse>
    Differentiable<ScalarType, DiffMode::Reverse>::random_normal(RandomGenerator& gen) {
        return Differentiable(ScalarType::random_normal(gen));
    }

    template<class ScalarType>
    template<class Distribution, class RandomGenerator>
    inline Differentiable<ScalarType, DiffMode::Reverse>
    Differentiable<ScalarType, DiffMode::Reverse>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Differentiable(ScalarType::random_any(dist, gen));
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator+(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() + s2.getValue(), s1.getTangent() + s2.getTangent());
            else
                return ResultType(s1.getValue() + s2.getValue(), s1.getTangent());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType>::getInstance();
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
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        return s2 + s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator-(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() - s2.getValue(), s1.getTangent() - s2.getTangent());
            else
                return ResultType(s1.getValue() - s2.getValue(), s1.getTangent());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType>::getInstance();
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
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        return -(s2 - s1.getDerived());
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator*(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() * s2.getValue(), s1.getTangent() * s2.getValue() + s1.getValue() * s2.getTangent());
            else
                return ResultType(s1.getValue() * s2.getValue(), s1.getTangent() * s2.getValue());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType>::getInstance();
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
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        return s2 * s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator/(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using FirstType = Differentiable<ScalarType, Mode>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            const auto rep = reciprocal(s2.getValue());
            if constexpr (OtherScalar::isDifferentiable)
                return ResultType(s1.getValue() * rep, (s1.getTangent() * s2.getValue() - s1.getValue() * s2.getTangent()) * square(rep));
            else
                return ResultType(s1.getValue() * rep, s1.getTangent() * rep);
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType>::getInstance();
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
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            const auto rep = reciprocal(s2.getValue());
            return ResultType(s1.getValue() * rep, -s1.getValue() * s2.getTangent() * square(rep));
        }
        else {
            static_assert(std::is_same<OtherScalar, ScalarType>::value, "[Error]: Reverse mode between different type is not supported");
            const ScalarType value = s1.getValue() / s2.getValue();
            const ResultType copy = s1.getDerived();
            auto& tracer = DiffTracer<ScalarType>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Div);
            tracer.pushOperand(copy);
            tracer.pushOperand(s2.getDerived());
            return result;
        }
    }
}
