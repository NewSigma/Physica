/*
 * Copyright 2024-2025 Weibo He.
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

#include "../VectorExpr.h"

namespace Physica {
    template<Vector V, Scalar U>
    class VectorExpr<ExprType::Mul, V, U>
            : public BinaryVectorExpr<ExprType::Mul, V, U> {
        using Base = BinaryVectorExpr<ExprType::Mul, V, U>;
    public:
        using Base::isReverseDiff;
        using Base::isComplex;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] auto operator-() const& noexcept;
        [[nodiscard]] auto operator-() && noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        inline void assign(Vector auto& v) const;

        template<ExecutePolicy P = Sequential>
        inline void assign_add(Vector auto& v) const;
        void assign_add_base(Vector auto& v) const noexcept;
        void assign_add_mkl(void* y) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        [[nodiscard]] T sum() const { return getLHS().sum() * getRHS(); }

        void reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff);
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V, Scalar U>
    auto VectorExpr<ExprType::Mul, V, U>::operator-() const& noexcept {
        return getLHS() * (-getRHS());
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprType::Mul, V, U>::operator-() && noexcept {
        return std::move(getLHS()) * (-getRHS());
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    inline void VectorExpr<ExprType::Mul, V, U>::assign(Vector auto& v) const {
        constexpr bool FastAssign = Traits<std::remove_cvref_t<V>>::FastAssign;
        if constexpr (FastAssign) {
            getLHS().template assign<P>(v);
            v *= getRHS();
        }
        else
            Base::assign(v);
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    inline void VectorExpr<ExprType::Mul, V, U>::assign_add(Vector auto& v) const {
        using V1 = std::remove_cvref_t<decltype(v)>;
        constexpr size_t Size = std::max(Base::SizeAtCompile, V1::SizeAtCompile);
        constexpr bool SmallVector = 0 < Size && Size <= 128;
        if constexpr (Internal::EnableMKL<V, V1>::value && !SmallVector) {
            if (Base::getLength() <= 128)
                assign_add_base(v);
            else
                assign_add_mkl(v.data());
        }
        else
            assign_add_base(v);
    }

    template<Vector V, Scalar U>
    void VectorExpr<ExprType::Mul, V, U>::assign_add_base(Vector auto& v) const noexcept {
        Base::assign_add(v);
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprType::Mul, V, U>::calc(size_t index) const -> CoDiff<T> {
        return getLHS().calc(index) * getRHS();
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprType::Mul, V, U>::calc_value(size_t index) const -> Tv {
        return getLHS().calc_value(index) * getRHS().value();
    }

    template<Vector V, Scalar U>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Mul, V, U>::packet(size_t index) const {
        return getLHS().template packet<Pack>(index) * Pack(getRHS());
    }

    template<Vector V, Scalar U>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Mul, V, U>::packetPartial(size_t index, size_t count) const {
        return getLHS().template packetPartial<Pack>(index, count) * Pack(getRHS());
    }

    template<Vector V, Scalar U>
    void VectorExpr<ExprType::Mul, V, U>::reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff) {
        const auto& grad = grad_.values();
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (ReverseDiff<V>)
            lhs.reverse(rhs.value() * grad);
        if constexpr (ReverseDiff<U>)
            rhs.reverse(lhs.values() * grad);
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprType::Mul, V1, V2>
            : public BinaryVectorExpr<ExprType::Mul, V1, V2> {
        using Base = BinaryVectorExpr<ExprType::Mul, V1, V2>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS().calc(index);
        }

        [[nodiscard]] Tv calc_value(size_t index) const {
            return Base::getLHS().calc_value(index) * Base::getRHS().calc_value(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index) * Base::getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count) * Base::getRHS().template packetPartial<Pack>(index, count);
        }

        void reverse(const auto& grad_) const noexcept requires(isReverseDiff);
    };

    template<Vector V1, Vector V2>
    void VectorExpr<ExprType::Mul, V1, V2>::reverse(const auto& grad_) const noexcept requires(isReverseDiff) {
        if constexpr (Scalar<decltype(grad_)>) {
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(Base::getRHS().values() * grad_);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(Base::getLHS().values() * grad_);
        }
        else {
            static_assert(Vector<decltype(grad_)>);
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(hadamard(Base::getRHS().values(), grad));
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(hadamard(Base::getLHS().values(), grad));
        }
    }

    template<Vector V, Scalar U>
    [[nodiscard]] inline auto operator*(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Mul, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] inline auto operator*(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return std::forward<V>(v) * std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] inline auto hadamard(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprType::Mul, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#include "MKL/Mul.h"
