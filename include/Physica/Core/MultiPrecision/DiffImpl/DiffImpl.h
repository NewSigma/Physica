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

namespace Physica::Core {
    template<class ScalarType, int Order>
    Diff<ScalarType, DiffMode::Forward, Order>::Diff(ScalarType value_)
            : value(std::move(value_)), grad(0) {}

    template<class ScalarType, int Order>
    Diff<ScalarType, DiffMode::Forward, Order>::Diff(ScalarType value_, GradType grad_)
            : value(std::move(value_)), grad(std::move(grad_)) {}

    template<class ScalarType, int Order>
    template<class OtherScalar, int OtherOrder>
    Diff<ScalarType, DiffMode::Forward, Order>::Diff(const Diff<OtherScalar, DiffMode::Forward, OtherOrder>& other)
            : value(other.getValue()), grad(other.getGrad()) {}

    template<class ScalarType, int Order>
    inline bool Diff<ScalarType, DiffMode::Forward, Order>::operator==(const This& other) const {
        return value == other.value && grad == other.grad;
    }

    template<class ScalarType, int Order>
    inline Diff<ScalarType, DiffMode::Forward, Order> Diff<ScalarType, DiffMode::Forward, Order>::operator-() const {
        return {-value, -grad};
    }

    template<class ScalarType, int Order>
    Diff<ScalarType, DiffMode::Forward, Order> Diff<ScalarType, DiffMode::Forward, Order>::conjugate() const {
        return This(getValue().conjugate(), getGrad().conjugate());
    }

    template<class ScalarType, int Order>
    void Diff<ScalarType, DiffMode::Forward, Order>::swap(Diff& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        value.swap(obj.value);
        grad.swap(obj.grad);
    }

    template<class ScalarType, int Order>
    template<int GradOrder>
    typename Diff<ScalarType, DiffMode::Forward, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<ScalarType, DiffMode::Forward, Order>::getGrad() noexcept {
        static_assert(Order >= GradOrder, "[Error]: Order is not enough to calculate the required grad");
        static_assert(GradOrder > 0, "[Error]: 0 or minus order is not well defined");
        if constexpr (GradOrder == 1)
            return grad;
        else
            return getGrad().template getGrad<GradOrder - 1>();
    }

    template<class ScalarType, int Order>
    template<int GradOrder>
    inline const typename Diff<ScalarType, DiffMode::Forward, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<ScalarType, DiffMode::Forward, Order>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<class ScalarType, int Order>
    template<class RandomGenerator>
    inline Diff<ScalarType, DiffMode::Forward, Order>
    Diff<ScalarType, DiffMode::Forward, Order>::random_uniform(RandomGenerator& gen) {
        return Diff(ScalarType::random_uniform(gen));
    }

    template<class ScalarType, int Order>
    template<class RandomGenerator>
    inline Diff<ScalarType, DiffMode::Forward, Order>
    Diff<ScalarType, DiffMode::Forward, Order>::random_normal(RandomGenerator& gen) {
        return Diff(ScalarType::random_normal(gen));
    }

    template<class ScalarType, int Order>
    template<class Distribution, class RandomGenerator>
    [[nodiscard]] inline Diff<ScalarType, DiffMode::Forward, Order>
    Diff<ScalarType, DiffMode::Forward, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Diff(ScalarType::random_any(dist, gen));
    }

#ifdef PHYSICA_HDF5
    template<class ScalarType, int Order>
    const H5::DataType& Diff<ScalarType, DiffMode::Forward, Order>::getH5DataType() {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Value", HOFFSET(This, value), ScalarType::getH5DataType().getId());
            H5Tinsert(id, "Grad", HOFFSET(This, grad), GradType::getH5DataType().getId());
            return result;
        }());
        return *instance;
    }
#endif

    template<class T, int Order>
    ScalarRef<Diff<T, DiffMode::Forward, Order>>& ScalarRef<Diff<T, DiffMode::Forward, Order>>::operator=(const ScalarRef& other) {
        getValue() = other.getValue();
        getGrad() = other.getGrad();
        return *this;
    }

    template<class T, int Order>
    ScalarRef<Diff<T, DiffMode::Forward, Order>>& ScalarRef<Diff<T, DiffMode::Forward, Order>>::operator=(const ScalarType& other) {
        getValue() = other.getValue();
        getGrad() = other.getGrad();
        return *this;
    }

    template<class T, int Order>
    void ScalarRef<Diff<T, DiffMode::Forward, Order>>::swap(This&& obj) noexcept {
        getValue().swap(obj.getValue());
        getGrad().swap(obj.getGrad());
    }

    template<class T, int Order>
    void ScalarRef<Diff<T, DiffMode::Forward, Order>>::swap(ScalarType& obj) noexcept {
        getValue().swap(obj.getValue());
        getGrad().swap(obj.getGrad());
    }

    template<class T, int Order>
    template<int GradOrder>
    typename ScalarRef<Diff<T, DiffMode::Forward, Order>>::template GradRtnTy<GradOrder> ScalarRef<Diff<T, DiffMode::Forward, Order>>::getGrad() noexcept {
        if constexpr (GradOrder == 1) {
            if constexpr (Order == GradOrder)
                return *ptr.grad_ptr();
            else
                return ScalarRef<Diff<T, DiffMode::Forward, Order - GradOrder>>(ptr.grad_ptr());
        }
        else
            return ptr.grad_ptr().template getGrad<GradOrder - 1>();
    }

    template<class T, int Order>
    template<int GradOrder>
    const typename ScalarRef<Diff<T, DiffMode::Forward, Order>>::template GradRtnTy<GradOrder> ScalarRef<Diff<T, DiffMode::Forward, Order>>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType, int Order>
    Diff<ScalarType, DiffMode::Reverse, Order>::Diff(ScalarType value)
            : Diff(TracerType::getInstance().pushOperation(std::move(value), ExprType::Set)) {}

    template<class ScalarType, int Order>
    Diff<ScalarType, DiffMode::Reverse, Order>::Diff(
            [[maybe_unused]] ScalarType value, [[maybe_unused]] ScalarType grad) {
        throw std::runtime_error("[Error]: This function is provided for template meta programming, you should not arrive here");
    }

    template<class ScalarType, int Order>
    Diff<ScalarType, DiffMode::Reverse, Order>::Diff(ValueType pValue_, GradType grad_)
            : pValue(pValue_), grad(grad_) {}

    template<class ScalarType, int Order>
    template<int KeepDeep>
    Diff<ScalarType, DiffMode::Reverse, Order>::operator Diff<ScalarType, DiffMode::Reverse, KeepDeep>() const {
        static_assert(KeepDeep < Order, "[Error]: Keep depth equal or greater than order makes no sense");
        static_assert(KeepDeep > 0, "[Error]: Keep 0 depth does nothing");
        using ResultType = Diff<ScalarType, DiffMode::Reverse, KeepDeep>;
        if constexpr (KeepDeep == 1)
            return ResultType(value_ptr(), grad_ptr());
        else
            return ResultType(value_ptr(), Diff<ScalarType, DiffMode::Reverse, KeepDeep - 1>(grad));
    }

    template<class ScalarType, int Order>
    inline Diff<ScalarType, DiffMode::Reverse, Order>& Diff<ScalarType, DiffMode::Reverse, Order>::operator=(std::nullptr_t) noexcept {
        pValue = nullptr;
        grad = nullptr;
        return *this;
    }

    template<class ScalarType, int Order>
    inline bool Diff<ScalarType, DiffMode::Reverse, Order>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (grad == other.grad) && "[Error]: Bad scalar");
        return result;
    }

    template<class ScalarType, int Order>
    inline bool Diff<ScalarType, DiffMode::Reverse, Order>::operator==(std::nullptr_t) const noexcept {
        const bool isInvalid = value_ptr() == nullptr;
        return isInvalid;
    }

    template<class ScalarType, int Order>
    inline Diff<ScalarType, DiffMode::Reverse, Order> Diff<ScalarType, DiffMode::Reverse, Order>::operator-() const {
        auto& tracer = TracerType::getInstance();
        const auto result = tracer.pushOperation(-getValue(), ExprType::Minus);
        tracer.pushOperand(*this);
        return result;
    }

    template<class ScalarType, int Order>
    inline void Diff<ScalarType, DiffMode::Reverse, Order>::reverse() {
        if (getSource() == ExprType::Diff) { //Optimize: Always false if no high order differential is required
            getFirstOperand()->reverse();
            return;
        }

        if constexpr (Order == 1)
            *grad = ScalarType(1);
        else
            grad.setValue(ScalarType(1));
        TracerType::getInstance().reverse_from(*this);
    }

    template<class ScalarType, int Order>
    inline void Diff<ScalarType, DiffMode::Reverse, Order>::reverse_to(Diff to) {
        if (getSource() == ExprType::Diff) {
            getFirstOperand()->reverse_to(to);
            return;
        }

        if constexpr (Order == 1)
            *grad = ScalarType(1);
        else
            grad.setValue(ScalarType(1));
        TracerType::getInstance().reverse(*this, to);
    }

    template<class ScalarType, int Order>
    inline void Diff<ScalarType, DiffMode::Reverse, Order>::zero_grad() {
        *grad_ptr() = ScalarType(0);
    }

    template<class ScalarType, int Order>
    inline Diff<ScalarType, DiffMode::Reverse, Order> Diff<ScalarType, DiffMode::Reverse, Order>::copy() const {
        return TracerType::getInstance().copy(*this);
    }

    template<class ScalarType, int Order>
    inline void Diff<ScalarType, DiffMode::Reverse, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(grad, obj.grad);
    }

    template<class ScalarType, int Order>
    __host__ __device__ inline ScalarType* Diff<ScalarType, DiffMode::Reverse, Order>::value_ptr() const noexcept {
        return pValue;
    }

    template<class ScalarType, int Order>
    __host__ __device__ inline ScalarType* Diff<ScalarType, DiffMode::Reverse, Order>::grad_ptr() const noexcept {
        if constexpr (Order == 1)
            return grad;
        else
            return grad.value_ptr();
    }

    template<class ScalarType, int Order>
    inline ScalarType& Diff<ScalarType, DiffMode::Reverse, Order>::getValue() noexcept {
        return *pValue;
    }

    template<class ScalarType, int Order>
    template<int GradOrder>
    inline typename Diff<ScalarType, DiffMode::Reverse, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<ScalarType, DiffMode::Reverse, Order>::getGrad() noexcept {
        static_assert(Order >= GradOrder, "[Error]: Order is not enough to calculate the required grad");
        static_assert(GradOrder > 0, "[Error]: 0 or minus order is not well defined");
        if constexpr (GradOrder == 1) {
            if constexpr (Order == 1)
                return *grad;
            else
                return grad;
        }
        else
            return getGrad<1>().template getGrad<GradOrder - 1>();
    }

    template<class ScalarType, int Order>
    inline const ScalarType& Diff<ScalarType, DiffMode::Reverse, Order>::getValue() const noexcept {
        return const_cast<This&>(*this).getValue();
    }

    template<class ScalarType, int Order>
    template<int GradOrder>
    inline const typename Diff<ScalarType, DiffMode::Reverse, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<ScalarType, DiffMode::Reverse, Order>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<class ScalarType, int Order>
    inline ExprType Diff<ScalarType, DiffMode::Reverse, Order>::getSource() const noexcept {
        auto& tracer = TracerType::getInstance();
        return tracer.getSource(*this);
    }

    template<class ScalarType, int Order>
    inline Diff<ScalarType, DiffMode::Reverse, Order>*
    Diff<ScalarType, DiffMode::Reverse, Order>::getFirstOperand() const {
        auto& tracer = TracerType::getInstance();
        return tracer.getFirstOperand(*this);
    }

    template<class ScalarType, int Order>
    template<class RandomGenerator>
    inline Diff<ScalarType, DiffMode::Reverse, Order>
    Diff<ScalarType, DiffMode::Reverse, Order>::random_uniform(RandomGenerator& gen) {
        return Diff(ScalarType::random_uniform(gen));
    }

    template<class ScalarType, int Order>
    template<class RandomGenerator>
    inline Diff<ScalarType, DiffMode::Reverse, Order>
    Diff<ScalarType, DiffMode::Reverse, Order>::random_normal(RandomGenerator& gen) {
        return Diff(ScalarType::random_normal(gen));
    }

    template<class ScalarType, int Order>
    template<class Distribution, class RandomGenerator>
    inline Diff<ScalarType, DiffMode::Reverse, Order>
    Diff<ScalarType, DiffMode::Reverse, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Diff(ScalarType::random_any(dist, gen));
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline auto operator+(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2_) {
        using FirstType = Diff<ScalarType, Mode, Order>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        const auto& s2 = s2_.getDerived();
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
            
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const PlainScalar value = s1.getValue() + s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExprType::Add);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExprType::Add);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2) {
        return s2 + s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline auto operator-(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2_) {
        using FirstType = Diff<ScalarType, Mode, Order>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        const auto& s2 = s2_.getDerived();
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
            
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const PlainScalar value = s1.getValue() - s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExprType::Sub);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExprType::Sub);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2) {
        return -(s2 - s1.getDerived());
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline auto operator*(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2_) {
        using FirstType = Diff<ScalarType, Mode, Order>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        const auto& s2 = s2_.getDerived();
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable) {
                using GradType = typename ResultType::GradType;
                return ResultType(s1.getValue() * s2.getValue(), GradType(GradType(s2) * s1.getGrad() + GradType(s1) * s2.getGrad()));
            }
            else
                return ResultType(s1.getValue() * s2.getValue(), s1.getGrad() * s2.getValue());
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const PlainScalar value = s1.getValue() * s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExprType::Mul);
                tracer.pushOperand(s1);
                tracer.pushOperand(s2);
            }
            else {
                const ResultType copy = s2;
                result = tracer.pushOperation(value, ExprType::Mul);
                tracer.pushOperand(s1);
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2) {
        return s2 * s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline auto operator/(const Diff<ScalarType, Mode, Order>& s1, const ScalarBase<OtherScalar>& s2_) {
        using FirstType = Diff<ScalarType, Mode, Order>;
        using SecondType = OtherScalar;
        using ResultType = typename Internal::BinaryScalarOpReturnType<FirstType, SecondType>::Type;
        const auto& s2 = s2_.getDerived();
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (OtherScalar::isDifferentiable) {
                using GradType = typename ResultType::GradType;
                const auto v = reciprocal(GradType(s2));
                return ResultType(s1.getValue() * v.getValue(), GradType((s1.getGrad() * GradType(s2) - GradType(s1) * s2.getGrad()) * square(v)));
            }
            else {
                const auto v = reciprocal(s2.getValue());
                return ResultType(s1.getValue() * v, s1.getGrad() * v);
            }
        }
        else {
            using FirstPlain = typename FirstType::PlainScalar;
            using SecondPlain = typename SecondType::PlainScalar;
            using PlainScalar = typename Internal::BinaryScalarOpReturnType<FirstPlain, SecondPlain>::Type;
            static_assert(std::is_same<FirstPlain, SecondPlain>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const PlainScalar value = s1.getValue() / s2.getValue();
            ResultType result;
            if constexpr (OtherScalar::isDifferentiable) {
                result = tracer.pushOperation(value, ExprType::Div);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(s2.getDerived());
            }
            else {
                const ResultType copy = s2.getDerived();
                result = tracer.pushOperation(value, ExprType::Div);
                tracer.pushOperand(s1.getDerived());
                tracer.pushOperand(copy);
            }
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Diff<ScalarType, Mode, Order>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Diff<ScalarType, Mode, Order>, OtherScalar>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            const auto rep = reciprocal(s2.getValue());
            return ResultType(s1.getValue() * rep, -s1.getValue() * s2.getGrad() * square(rep));
        }
        else {
            static_assert(std::is_same<OtherScalar, ScalarType>::value, "[Error]: Reverse mode between different type is not supported");
            const ScalarType value = s1.getValue() / s2.getValue();
            const ResultType copy = s1.getDerived();
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Div);
            tracer.pushOperand(copy);
            tracer.pushOperand(s2.getDerived());
            return result;
        }
    }
}
