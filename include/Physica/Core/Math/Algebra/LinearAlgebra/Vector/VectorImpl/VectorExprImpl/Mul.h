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
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto&& v) const;
        void assign_add_mkl(Vector auto&& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add_base(Vector auto&& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        [[nodiscard]] T sum() const { return getLHS().sum() * getRHS(); }

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
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
    void VectorExpr<ExprType::Mul, V, U>::assign(Vector auto&& v) const {
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
    void VectorExpr<ExprType::Mul, V, U>::assign_add(Vector auto&& v) const {
        using V1 = std::remove_cvref_t<decltype(v)>;
        constexpr size_t Length = std::max(Base::SizeAtCompile, V1::SizeAtCompile);
        constexpr bool SmallVector = 0 < Length && Length <= 128;
        if constexpr (Internal::EnableMKL<V, V1>::value && !SmallVector) {
            if (Base::getLength() <= 128)
                assign_add_base(v);
            else
                assign_add_mkl(v);
        }
        else
            assign_add_base(v);
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Mul, V, U>::assign_add_base(Vector auto&& v) const noexcept {
        Base::template assign_add<P>(v);
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
    void VectorExpr<ExprType::Mul, V, U>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (ReverseDiff<V>)
            lhs.reverse(rhs.value() * g);
        if constexpr (ReverseDiff<U>)
            rhs.reverse(lhs.values() * g);
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprType::Mul, V1, V2>
            : public BinaryVectorExpr<ExprType::Mul, V1, V2> {
        using Base = BinaryVectorExpr<ExprType::Mul, V1, V2>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        using Base::reverse;
        void reverse(const auto& grad) const noexcept;
    };

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Mul, V1, V2>::assign(Vector auto&& v) const {
        Base::template assign_base<P>(v);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprType::Mul, V1, V2>::calc(size_t index) const -> CoDiff<T> {
        return Base::getLHS().calc(index) * Base::getRHS().calc(index);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprType::Mul, V1, V2>::calc_value(size_t index) const -> Tv {
        return Base::getLHS().calc_value(index) * Base::getRHS().calc_value(index);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Mul, V1, V2>::packet(size_t index) const {
        return Base::getLHS().template packet<Pack>(index) * Base::getRHS().template packet<Pack>(index);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Mul, V1, V2>::packetPartial(size_t index, size_t count) const {
        return Base::getLHS().template packetPartial<Pack>(index, count) * Base::getRHS().template packetPartial<Pack>(index, count);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprType::Mul, V1, V2>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        if constexpr (Scalar<decltype(grad)>) {
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(Base::getRHS().values() * grad);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(Base::getLHS().values() * grad);
        }
        else {
            static_assert(Vector<decltype(grad)>);
            const auto& g = grad.values();
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(hadamard(Base::getRHS().values(), g));
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(hadamard(Base::getLHS().values(), g));
        }
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator*(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        using RtnTy = VectorExpr<ExprType::Mul, V&&, U&&>;
        if constexpr (instanceof_xt<VectorExpr, V>) {
            constexpr ExprType ID = Traits<V>::Type;
            using RHS = Traits<V>::RHS;
            if constexpr (ID == ExprType::Mul && Scalar<RHS>)
                return v.getLHS() * (v.getRHS() * x);
            else
                return RtnTy(std::forward<V>(v), std::forward<U>(x));
        }
        else
            return RtnTy(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] auto operator*(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return std::forward<V>(v) * std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto hadamard(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        using RtnTy = VectorExpr<ExprType::Mul, V1&&, V2&&>;
        if constexpr (instanceof_xt<VectorExpr, V1>) {
            using RHS1 = Traits<V1>::RHS;
            constexpr ExprType ID1 = Traits<V1>::Type;
            if constexpr (ID1 == ExprType::Mul && Scalar<RHS1>) // if we can lower V1
                return hadamard(v2, v1.getLHS()) * v1.getRHS();
            else if constexpr (instanceof_xt<VectorExpr, V2>) { // if not, see if we can lower V2
                using RHS2 = Traits<V2>::RHS;
                constexpr ExprType ID2 = Traits<V2>::Type;
                if constexpr (ID2 == ExprType::Mul && Scalar<RHS2>)
                    return hadamard(v1, v2.getLHS()) * v2.getRHS();
                else
                    return RtnTy(std::forward<V1>(v1), std::forward<V2>(v2));
            }
            else // there is nothing we can do here
                return RtnTy(std::forward<V1>(v1), std::forward<V2>(v2));
        }
        else if constexpr (instanceof_xt<VectorExpr, V2>)
            return hadamard(std::forward<V2>(v2), std::forward<V1>(v1));
        else
            return RtnTy(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#include "MKL/Mul.h"
