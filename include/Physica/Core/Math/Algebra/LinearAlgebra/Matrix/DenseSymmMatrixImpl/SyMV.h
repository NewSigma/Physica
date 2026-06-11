/*
 * Copyright 2020-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using typename Base::Tv;
        using typename Base::Trv;
    private:
        decay_rvalue_t<M> m;
        decay_rvalue_t<V> v;
    public:
        GEMV(M&& m_, V&& v_);
        GEMV(const This&) = delete;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return m.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    GEMV<M, V>::GEMV(M&& m_, V&& v_) : m(m_), v(v_) {
        assert(m.getCol() == v.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        const size_t length = getLength();
        assert(length == target.getLength());
        if (length >= 16) {
            target.zeros();
            for (size_t i = 0; i < length; ++i) {
                const size_t diag = m.toIndex1D(i, i);
                const auto seg1 = m.asVector().segment(diag, diag + length - i);
                const auto seg2 = v.segment(i, length);
                target[i] += seg1 * seg2;

                if (i + 1 < length)
                    (seg1.tail(1) * v.calc(i)).assign_add(target.segment(i + 1, length));
            }
        }
        else
            Base::assign(target);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    auto GEMV<M, V>::calc(size_t index) const -> CoDiff<ScalarType> {
        return m.row(index) * v;
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    auto GEMV<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.m);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, DenseSymmMatrix>)
    __host__ __device__ consteval size_t GEMV<M, V>::getSizeAtCompile() noexcept {
        using M1 = std::remove_cvref<M>::type;
        using V1 = std::remove_cvref<V>::type;
        return std::max(M1::getRowAtCompile(), V1::getSizeAtCompile());
    }
}
