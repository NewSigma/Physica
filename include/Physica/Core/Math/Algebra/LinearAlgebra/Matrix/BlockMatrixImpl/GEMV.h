/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/BlockMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<BlockMatrix, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    private:
        LazyDestroy<M> m;
        LazyDestroy<V> v;
    public:
        GEMV(M m_, V v_);
        GEMV(const This&) = delete;
        GEMV(This&&) noexcept = delete;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target_) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        [[nodiscard]] Tv calc_value(size_t) const { noImpl(__func__); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof<BlockMatrix, M>)
    GEMV<M, V>::GEMV(M m_, V v_) : m(std::forward<M>(m_)), v(std::forward<V>(v_)) {
        assert(m.getCol() == v.getLength() && "[Error]: Dimensions do not match");
    }

    template<Matrix M, Vector V> requires(instanceof<BlockMatrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        assert(getLength() == target.getLength() && "[Error]: Dimensions do not match");
        size_t from = 0;
        for (size_t i = 0; i < m.getNumBlocks(); ++i) {
            const size_t to = m.getIndexEnds()[i];
            const auto v1 = v.segment(from, to);
            auto target1 = target.segment(from, to);
            target1 = m.getBlocks()[i] * v1;
            from = to;
        }
    }

    template<Matrix M, Vector V> requires(instanceof<BlockMatrix, M>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.m);
    }

    template<Matrix M, Vector V> requires(instanceof<BlockMatrix, M>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v);
    }
}
