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
    template<Scalar T, int Order>
    Diff<T, DiffMode::Forward, Order>::Diff(T value_)
            : value(std::move(value_)), grad(0) {}

    template<Scalar T, int Order>
    Diff<T, DiffMode::Forward, Order>::Diff(T value_, GradType grad_)
            : value(std::move(value_)), grad(std::move(grad_)) {}

    template<Scalar T, int Order>
    template<Scalar U, int OtherOrder>
    Diff<T, DiffMode::Forward, Order>::Diff(const Diff<U, DiffMode::Forward, OtherOrder>& other)
            : value(other.getValue()), grad(other.getGrad()) {}

    template<Scalar T, int Order>
    inline bool Diff<T, DiffMode::Forward, Order>::operator==(const This& other) const {
        return value == other.value && grad == other.grad;
    }

    template<Scalar T, int Order>
    inline Diff<T, DiffMode::Forward, Order> Diff<T, DiffMode::Forward, Order>::operator-() const {
        return {-value, -grad};
    }

    template<Scalar T, int Order>
    Diff<T, DiffMode::Forward, Order> Diff<T, DiffMode::Forward, Order>::conjugate() const {
        return This(getValue().conjugate(), getGrad().conjugate());
    }

    template<Scalar T, int Order>
    void Diff<T, DiffMode::Forward, Order>::swap(Diff& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        value.swap(obj.value);
        grad.swap(obj.grad);
    }

    template<Scalar T, int Order>
    void Diff<T, DiffMode::Forward, Order>::swap(ScalarRef<This>&& ref) noexcept {
        assert(ScalarPtr<This>(*this) != ScalarPtr<This>(ref) && "[Error]: Self swap is likely a bug");
        value.swap(ref.getValue());
        grad.swap(ref.getGrad());
    }

    template<Scalar T, int Order>
    template<int GradOrder>
    Diff<T, DiffMode::Forward, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<T, DiffMode::Forward, Order>::getGrad() noexcept {
        static_assert(Order >= GradOrder, "[Error]: Order is not enough to calculate the required grad");
        static_assert(GradOrder > 0, "[Error]: 0 or minus order is not well defined");
        if constexpr (GradOrder == 1)
            return grad;
        else
            return getGrad().template getGrad<GradOrder - 1>();
    }

    template<Scalar T, int Order>
    template<int GradOrder>
    inline const Diff<T, DiffMode::Forward, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<T, DiffMode::Forward, Order>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<Scalar T, int Order>
    __host__ __device__ inline bool Diff<T, DiffMode::Forward, Order>::isFinite() const noexcept {
        return value.isFinite() && grad.isFinite();
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline Diff<T, DiffMode::Forward, Order>
    Diff<T, DiffMode::Forward, Order>::random_uniform() {
        return Diff(T::template random_uniform<R>());
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline Diff<T, DiffMode::Forward, Order>
    Diff<T, DiffMode::Forward, Order>::random_normal() {
        return Diff(T::template random_normal<R>());
    }

    template<Scalar T, int Order>
    template<class Distribution, RandomGenerator R>
    [[nodiscard]] inline Diff<T, DiffMode::Forward, Order>
    Diff<T, DiffMode::Forward, Order>::random_any(Distribution& dist) {
        return Diff(T::template random_any<R>(dist));
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, int Order>
    const H5::DataType& Diff<T, DiffMode::Forward, Order>::getH5DataType() {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Value", HOFFSET(This, value), T::getH5DataType().getId());
            H5Tinsert(id, "Grad", HOFFSET(This, grad), GradType::getH5DataType().getId());
            return result;
        }());
        return *instance;
    }
#endif

    template<Scalar T, int Order>
    inline ScalarRef<Diff<T, DiffMode::Forward, Order>>& ScalarRef<Diff<T, DiffMode::Forward, Order>>::operator=(const ScalarRef& other) {
        getValue() = other.getValue();
        getGrad() = other.getGrad();
        return *this;
    }

    template<Scalar T, int Order>
    inline ScalarRef<Diff<T, DiffMode::Forward, Order>>& ScalarRef<Diff<T, DiffMode::Forward, Order>>::operator=(const ScalarType& other) {
        getValue() = other.getValue();
        getGrad() = other.getGrad();
        return *this;
    }

    template<Scalar T, int Order>
    void ScalarRef<Diff<T, DiffMode::Forward, Order>>::swap(This&& obj) noexcept {
        getValue().swap(obj.getValue());
        getGrad().swap(obj.getGrad());
    }

    template<Scalar T, int Order>
    void ScalarRef<Diff<T, DiffMode::Forward, Order>>::swap(ScalarType& obj) noexcept {
        obj.swap(*this);
    }

    template<Scalar T, int Order>
    template<int GradOrder>
    ScalarRef<Diff<T, DiffMode::Forward, Order>>::template GradRtnTy<GradOrder> ScalarRef<Diff<T, DiffMode::Forward, Order>>::getGrad() noexcept {
        if constexpr (GradOrder == 1) {
            if constexpr (Order == GradOrder)
                return *ptr.grad_ptr();
            else
                return ScalarRef<Diff<T, DiffMode::Forward, Order - GradOrder>>(ptr.grad_ptr());
        }
        else
            return (*ptr.grad_ptr()).template getGrad<GradOrder - 1>();
    }

    template<Scalar T, int Order>
    template<int GradOrder>
    const ScalarRef<Diff<T, DiffMode::Forward, Order>>::template GradRtnTy<GradOrder> ScalarRef<Diff<T, DiffMode::Forward, Order>>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<Scalar T, int Order>
    ScalarPtr<Diff<T, DiffMode::Forward, Order>>::ScalarPtr(ScalarType& x) : ScalarPtr(&x.getValue(), &x.getGrad()) {}

    template<Scalar T, int Order>
    ScalarPtr<Diff<T, DiffMode::Forward, Order>>::ScalarPtr(ScalarRef<ScalarType>& x) : ScalarPtr(&x.getValue(), &x.getGrad()) {}

    template<Scalar T, int Order>
    inline bool ScalarPtr<Diff<T, DiffMode::Forward, Order>>::operator==(const This& other) const noexcept {
        bool flag = pair.first == other.pair.first;
        assert(flag == (pair.second == other.pair.second) && "[Error]: Bad ScalarPtr");
        return flag;
    }

    template<Scalar T, int Order>
    inline ScalarPtr<Diff<T, DiffMode::Forward, Order>>& ScalarPtr<Diff<T, DiffMode::Forward, Order>>::operator++() {
        for (auto& p : arr)
            p++;
        return *this;
    }

    template<Scalar T, int Order>
    inline ScalarPtr<Diff<T, DiffMode::Forward, Order>>& ScalarPtr<Diff<T, DiffMode::Forward, Order>>::operator--() {
        for (auto& p : arr)
            p--;
        return *this;
    }

    template<Scalar T, int Order>
    inline const ScalarPtr<Diff<T, DiffMode::Forward, Order>> ScalarPtr<Diff<T, DiffMode::Forward, Order>>::operator++(int) {
        return This(pair.first++, pair.second++);
    }

    template<Scalar T, int Order>
    inline const ScalarPtr<Diff<T, DiffMode::Forward, Order>> ScalarPtr<Diff<T, DiffMode::Forward, Order>>::operator--(int) {
        return This(pair.first--, pair.second--);
    }

    template<Scalar T, int Order>
    inline void ScalarPtr<Diff<T, DiffMode::Forward, Order>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pair.swap(obj.pair);
    }
    ////////////////////////////////////////////////////////////
    template<Scalar T, int Order>
    Diff<T, DiffMode::Reverse, Order>::Diff(T value)
            : Diff(TracerType::getInstance().pushOperation(std::move(value), ExprType::Set)) {}

    template<Scalar T, int Order>
    Diff<T, DiffMode::Reverse, Order>::Diff(
            [[maybe_unused]] T value, [[maybe_unused]] T grad) {
        throw std::runtime_error("[Error]: This function is provided for template meta programming, you should not arrive here");
    }

    template<Scalar T, int Order>
    Diff<T, DiffMode::Reverse, Order>::Diff(ValuePtr pValue_, GradType grad_)
            : pValue(pValue_), grad(grad_) {}

    template<Scalar T, int Order>
    template<int KeepDeep>
    Diff<T, DiffMode::Reverse, Order>::operator Diff<T, DiffMode::Reverse, KeepDeep>() const {
        static_assert(KeepDeep < Order, "[Error]: Keep depth equal or greater than order makes no sense");
        static_assert(KeepDeep > 0, "[Error]: Keep 0 depth does nothing");
        using ResultType = Diff<T, DiffMode::Reverse, KeepDeep>;
        if constexpr (KeepDeep == 1)
            return ResultType(value_ptr(), grad_ptr());
        else
            return ResultType(value_ptr(), Diff<T, DiffMode::Reverse, KeepDeep - 1>(grad));
    }

    template<Scalar T, int Order>
    inline Diff<T, DiffMode::Reverse, Order>& Diff<T, DiffMode::Reverse, Order>::operator=(std::nullptr_t) noexcept {
        pValue = nullptr;
        grad = nullptr;
        return *this;
    }

    template<Scalar T, int Order>
    inline bool Diff<T, DiffMode::Reverse, Order>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (grad == other.grad) && "[Error]: Bad scalar");
        return result;
    }

    template<Scalar T, int Order>
    inline bool Diff<T, DiffMode::Reverse, Order>::operator==(std::nullptr_t) const noexcept {
        const bool isInvalid = value_ptr() == nullptr;
        return isInvalid;
    }

    template<Scalar T, int Order>
    inline Diff<T, DiffMode::Reverse, Order> Diff<T, DiffMode::Reverse, Order>::operator-() const {
        auto& tracer = TracerType::getInstance();
        const auto result = tracer.pushOperation(-getValue(), ExprType::Minus);
        tracer.pushOperand(*this);
        return result;
    }

    template<Scalar T, int Order>
    inline void Diff<T, DiffMode::Reverse, Order>::reverse() {
        if (getSource() == ExprType::Diff) { //Optimize: Always false if no high order differential is required
            getFirstOperand()->reverse();
            return;
        }

        if constexpr (Order == 1)
            *grad = T(1);
        else
            grad.setValue(T(1));
        TracerType::getInstance().reverse_from(*this);
    }

    template<Scalar T, int Order>
    inline void Diff<T, DiffMode::Reverse, Order>::reverse_to(Diff to) {
        if (getSource() == ExprType::Diff) {
            getFirstOperand()->reverse_to(to);
            return;
        }

        if constexpr (Order == 1)
            *grad = T(1);
        else
            grad.setValue(T(1));
        TracerType::getInstance().reverse(*this, to);
    }

    template<Scalar T, int Order>
    inline void Diff<T, DiffMode::Reverse, Order>::zero_grad() {
        *grad_ptr() = T(0);
    }

    template<Scalar T, int Order>
    inline Diff<T, DiffMode::Reverse, Order> Diff<T, DiffMode::Reverse, Order>::copy() const {
        return TracerType::getInstance().copy(*this);
    }

    template<Scalar T, int Order>
    inline void Diff<T, DiffMode::Reverse, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(grad, obj.grad);
    }

    template<Scalar T, int Order>
    __host__ __device__ inline T* Diff<T, DiffMode::Reverse, Order>::value_ptr() const noexcept {
        return pValue;
    }

    template<Scalar T, int Order>
    __host__ __device__ inline T* Diff<T, DiffMode::Reverse, Order>::grad_ptr() const noexcept {
        if constexpr (Order == 1)
            return grad;
        else
            return grad.value_ptr();
    }

    template<Scalar T, int Order>
    inline T& Diff<T, DiffMode::Reverse, Order>::getValue() noexcept {
        return *pValue;
    }

    template<Scalar T, int Order>
    template<int GradOrder>
    inline Diff<T, DiffMode::Reverse, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<T, DiffMode::Reverse, Order>::getGrad() noexcept {
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

    template<Scalar T, int Order>
    inline const T& Diff<T, DiffMode::Reverse, Order>::getValue() const noexcept {
        return const_cast<This&>(*this).getValue();
    }

    template<Scalar T, int Order>
    template<int GradOrder>
    inline const Diff<T, DiffMode::Reverse, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<T, DiffMode::Reverse, Order>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<Scalar T, int Order>
    inline ExprType Diff<T, DiffMode::Reverse, Order>::getSource() const noexcept {
        auto& tracer = TracerType::getInstance();
        return tracer.getSource(*this);
    }

    template<Scalar T, int Order>
    inline Diff<T, DiffMode::Reverse, Order>*
    Diff<T, DiffMode::Reverse, Order>::getFirstOperand() const {
        auto& tracer = TracerType::getInstance();
        return tracer.getFirstOperand(*this);
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline Diff<T, DiffMode::Reverse, Order>
    Diff<T, DiffMode::Reverse, Order>::random_uniform() {
        return Diff(T::template random_uniform<R>());
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline Diff<T, DiffMode::Reverse, Order>
    Diff<T, DiffMode::Reverse, Order>::random_normal() {
        return Diff(T::template random_normal<R>());
    }

    template<Scalar T, int Order>
    template<class Distribution, RandomGenerator R>
    inline Diff<T, DiffMode::Reverse, Order>
    Diff<T, DiffMode::Reverse, Order>::random_any(Distribution& dist) {
        return Diff(T::template random_any<decltype(dist), R>(dist));
    }
    ////////////////////////////////////////////////////////////
    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline auto operator+(const Diff<T, Mode, Order>& s1, const U& s2) {
        using FirstType = Diff<T, Mode, Order>;
        using SecondType = U;
        using ResultType = Internal::BinaryScalarOpRtnTy<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (U::isDifferentiable)
                return ResultType(s1.getValue() + s2.getValue(), s1.getGrad() + s2.getGrad());
            else
                return ResultType(s1.getValue() + s2.getValue(), s1.getGrad());
        }
        else {
            using FirstValue = FirstType::ValueType;
            using SecondValue = SecondType::ValueType;
            using ValueType = Internal::BinaryScalarOpRtnTy<FirstValue, SecondValue>::Type;
            static_assert(std::is_same<FirstValue, SecondValue>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<T, Order>::getInstance();
            const ValueType value = s1.getValue() + s2.getValue();
            ResultType result;
            if constexpr (U::isDifferentiable) {
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

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline std::enable_if<!U::isDifferentiable, typename Internal::BinaryScalarOpRtnTy<Diff<T, Mode, Order>, U>::Type>::type
    operator+(const U& s1, const Diff<T, Mode, Order>& s2) {
        return s2 + s1.getDerived();
    }

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline auto operator-(const Diff<T, Mode, Order>& s1, const U& s2) {
        using FirstType = Diff<T, Mode, Order>;
        using SecondType = U;
        using ResultType = Internal::BinaryScalarOpRtnTy<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (U::isDifferentiable)
                return ResultType(s1.getValue() - s2.getValue(), s1.getGrad() - s2.getGrad());
            else
                return ResultType(s1.getValue() - s2.getValue(), s1.getGrad());
        }
        else {
            using FirstValue = FirstType::ValueType;
            using SecondValue = SecondType::ValueType;
            using ValueType = Internal::BinaryScalarOpRtnTy<FirstValue, SecondValue>::Type;
            static_assert(std::is_same<FirstValue, SecondValue>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<T, Order>::getInstance();
            const ValueType value = s1.getValue() - s2.getValue();
            ResultType result;
            if constexpr (U::isDifferentiable) {
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

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline std::enable_if<!U::isDifferentiable, typename Internal::BinaryScalarOpRtnTy<Diff<T, Mode, Order>, U>::Type>::type
    operator-(const U& s1, const Diff<T, Mode, Order>& s2) {
        return -(s2 - s1.getDerived());
    }

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline auto operator*(const Diff<T, Mode, Order>& s1, const U& s2) {
        using FirstType = Diff<T, Mode, Order>;
        using SecondType = U;
        using ResultType = Internal::BinaryScalarOpRtnTy<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (U::isDifferentiable) {
                using GradType = ResultType::GradType;
                return ResultType(s1.getValue() * s2.getValue(), GradType(GradType(s2) * s1.getGrad() + GradType(s1) * s2.getGrad()));
            }
            else
                return ResultType(s1.getValue() * s2.getValue(), s1.getGrad() * s2.getValue());
        }
        else {
            using FirstValue = FirstType::ValueType;
            using SecondValue = SecondType::ValueType;
            using ValueType = Internal::BinaryScalarOpRtnTy<FirstValue, SecondValue>::Type;
            static_assert(std::is_same<FirstValue, SecondValue>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<T, Order>::getInstance();
            const ValueType value = s1.getValue() * s2.getValue();
            ResultType result;
            if constexpr (U::isDifferentiable) {
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

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline std::enable_if<!U::isDifferentiable, typename Internal::BinaryScalarOpRtnTy<Diff<T, Mode, Order>, U>::Type>::type
    operator*(const U& s1, const Diff<T, Mode, Order>& s2) {
        return s2 * s1.getDerived();
    }

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline auto operator/(const Diff<T, Mode, Order>& s1, const U& s2) {
        using FirstType = Diff<T, Mode, Order>;
        using SecondType = U;
        using ResultType = Internal::BinaryScalarOpRtnTy<FirstType, SecondType>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            if constexpr (U::isDifferentiable) {
                using GradType = ResultType::GradType;
                const auto v = reciprocal(GradType(s2));
                return ResultType(s1.getValue() * v.getValue(), GradType((s1.getGrad() * GradType(s2) - GradType(s1) * s2.getGrad()) * square(v)));
            }
            else {
                const auto v = reciprocal(s2.getValue());
                return ResultType(s1.getValue() * v, s1.getGrad() * v);
            }
        }
        else {
            using FirstValue = FirstType::ValueType;
            using SecondValue = SecondType::ValueType;
            using ValueType = Internal::BinaryScalarOpRtnTy<FirstValue, SecondValue>::Type;
            static_assert(std::is_same<FirstValue, SecondValue>::value, "[Error]: Reverse mode between different type is not supported");
            
            auto& tracer = DiffTracer<T, Order>::getInstance();
            const ValueType value = s1.getValue() / s2.getValue();
            ResultType result;
            if constexpr (U::isDifferentiable) {
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

    template<Scalar T, DiffMode Mode, int Order, Scalar U>
    [[nodiscard]] inline std::enable_if<!U::isDifferentiable, typename Internal::BinaryScalarOpRtnTy<Diff<T, Mode, Order>, U>::Type>::type
    operator/(const U& s1, const Diff<T, Mode, Order>& s2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Diff<T, Mode, Order>, U>::Type;
        if constexpr (Mode == DiffMode::Forward) {
            const auto rep = reciprocal(s2.getValue());
            return ResultType(s1.getValue() * rep, -s1.getValue() * s2.getGrad() * square(rep));
        }
        else {
            static_assert(std::is_same<T, U>::value, "[Error]: Reverse mode between different type is not supported");
            const T value = s1.getValue() / s2.getValue();
            const ResultType copy = s1.getDerived();
            auto& tracer = DiffTracer<T, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Div);
            tracer.pushOperand(copy);
            tracer.pushOperand(s2.getDerived());
            return result;
        }
    }
}
