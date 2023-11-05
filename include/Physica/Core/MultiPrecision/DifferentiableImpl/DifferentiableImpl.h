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
    template<class AnyScalar>
    Differentiable<ScalarType, DiffMode::Forward>::Differentiable(const ScalarBase<AnyScalar>& s)
            : value(s.getValue()), tangent(0) {}

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
    size_t Differentiable<ScalarType, DiffMode::Forward>::getTraceIndex() const {
        throw std::runtime_error("[Error]: This function is provided for template meta programming, you should not arrive here");
    }

    template<class ScalarType>
    template<class Distribution, class RandomGenerator>
    [[nodiscard]] inline Differentiable<ScalarType, DiffMode::Forward>
    Differentiable<ScalarType, DiffMode::Forward>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Differentiable(ScalarType::random_any(dist, gen));
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType>
    template<class AnyScalar>
    Differentiable<ScalarType, DiffMode::Reverse>::Differentiable(const ScalarBase<AnyScalar>& s)
            : value(s.getValue()) {
        index = DiffTracerType::getInstance().pushOperation(value, ExpressionType::Set);
    }

    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Reverse>::Differentiable(ScalarType value_, size_t index_)
            : value(value_), index(index_) {
        assert(index <= DiffTracerType::getInstance().getRecords().getLength());
    }

    template<class ScalarType>
    Differentiable<ScalarType, DiffMode::Reverse>::Differentiable(ScalarType value_, ScalarType tangent)
            : value(std::move(value_)) {
        index = DiffTracerType::getInstance().pushOperation(value, tangent, ExpressionType::Set);
        throw std::runtime_error("[Error]: This function is provided for template meta programming, you should not arrive here");
    }

    template<class ScalarType>
    inline bool Differentiable<ScalarType, DiffMode::Reverse>::operator==(const This& other) const {
        return index == other.index;
    }

    template<class ScalarType>
    inline Differentiable<ScalarType, DiffMode::Reverse> Differentiable<ScalarType, DiffMode::Reverse>::operator-() const {
        auto& tracer = DiffTracerType::getInstance();
        ScalarType v = -value;
        const size_t i = tracer.pushOperation(v, ExpressionType::Minus);
        tracer.pushOperand(index);
        return Differentiable(v, i);
    }

    template<class ScalarType>
    inline void Differentiable<ScalarType, DiffMode::Reverse>::reverse() {
        DiffTracerType::getInstance().reverse(index);
    }

    template<class ScalarType>
    void Differentiable<ScalarType, DiffMode::Reverse>::swap(Differentiable<ScalarType, DiffMode::Reverse>& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        value.swap(obj.value);
        std::swap(index, obj.index);
    }

    template<class ScalarType>
    inline const ScalarType& Differentiable<ScalarType, DiffMode::Reverse>::getTangent() const noexcept {
        return DiffTracerType::getInstance()[index].tangent;
    }

    template<class ScalarType>
    inline void Differentiable<ScalarType, DiffMode::Reverse>::setValue(const ScalarType& x) {
        value = x;
        DiffTracerType::getInstance()[index].value = x;
    }

    template<class ScalarType>
    template<class Distribution, class RandomGenerator>
    [[nodiscard]] inline Differentiable<ScalarType, DiffMode::Reverse>
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
            size_t index;
            if constexpr (OtherScalar::isDifferentiable) {
                index = tracer.pushOperation(value, ExpressionType::Add);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(s2.getTraceIndex());
            }
            else {
                const ResultType copy = s2;
                index = tracer.pushOperation(value, ExpressionType::Add);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(copy.getTraceIndex());
            }
            return ResultType(value, index);
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
            size_t index;
            if constexpr (OtherScalar::isDifferentiable) {
                index = tracer.pushOperation(value, ExpressionType::Sub);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(s2.getTraceIndex());
            }
            else {
                const ResultType copy = s2;
                index = tracer.pushOperation(value, ExpressionType::Sub);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(copy.getTraceIndex());
            }
            return ResultType(value, index);
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
            size_t index;
            if constexpr (OtherScalar::isDifferentiable) {
                index = tracer.pushOperation(value, ExpressionType::Mul);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(s2.getTraceIndex());
            }
            else {
                const ResultType copy = s2;
                index = tracer.pushOperation(value, ExpressionType::Mul);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(copy.getTraceIndex());
            }
            return ResultType(value, index);
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
            size_t index;
            if constexpr (OtherScalar::isDifferentiable) {
                index = tracer.pushOperation(value, ExpressionType::Div);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(s2.getTraceIndex());
            }
            else {
                const ResultType copy = s2;
                index = tracer.pushOperation(value, ExpressionType::Div);
                tracer.pushOperand(s1.getTraceIndex()).pushOperand(copy.getTraceIndex());
            }
            return ResultType(value, index);
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
            const ResultType copy = s1;
            auto& tracer = DiffTracer<ScalarType>::getInstance();
            const size_t index = tracer.pushOperation(value, ExpressionType::Div);
            tracer.pushOperand(copy.getTraceIndex()).pushOperand(s2.getTraceIndex());
            return ResultType(value, index);
        }
    }
}
