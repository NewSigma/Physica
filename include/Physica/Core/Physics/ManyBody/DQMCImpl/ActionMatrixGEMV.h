/*
 * Copyright 2026 Weibo He.
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

#include "ActionMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M&& mat, V&& vec);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto& target) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const { noImpl(__func__); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    private:
        void assign_add_potential(Vector auto& target) const noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    GEMV<M, V>::GEMV(M&& mat, V&& vec) : mat(std::forward<M>(mat)), vec(std::forward<V>(vec)) {}

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const noexcept {
        target.assert_assign(*this);

        const int numSite = mat.getNumSite();
        const int numFreq2 = mat.getNumFreq() * 2;
        (kronecker(IdentityMatrix<Trv>(numSite) * T(0, 1), mat.matsubara) * getRHS()).assign(target);
        (kronecker(mat.params.getHoppingMatrix() * mat.getBeta(), IdentityMatrix<Trv>(numFreq2)) * getRHS()).assign_add(target);
        assign_add_potential(target);
    }

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_add(Vector auto& target) const noexcept {
        target.assert_assign(*this);

        const int numSite = mat.getNumSite();
        const int numFreq2 = mat.getNumFreq() * 2;
        (kronecker(IdentityMatrix<Trv>(numSite) * T(0, 1), mat.matsubara) * getRHS()).assign_add(target);
        (kronecker(mat.params.getHoppingMatrix() * mat.getBeta(), IdentityMatrix<Trv>(numFreq2)) * getRHS()).assign_add(target);
        assign_add_potential(target);
    }

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    auto GEMV<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    void GEMV<M, V>::assign_add_potential(Vector auto& target) const noexcept {
        const int numSite = mat.getNumSite();
        const int numFreq2 = mat.getNumFreq() * 2;
        auto buffer = MatrixND<T>(numFreq2);
        buffer.zeros();
        for (int site = 0, from = 0; site < numSite; ++site) {
            int to = from + numFreq2;
            mat.assign_potential(buffer, site);
            (buffer * getRHS().segment(from, to)).assign_add(target.segment(from, to));
            from += numFreq2;
        }
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<ActionMatrix, M>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        static_assert(M1::ColAtCompile == V1::SizeAtCompile || M1::ColAtCompile == Dynamic || V1::SizeAtCompile == Dynamic,
                "Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename V1::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = M1::RowAtCompile;
        constexpr static bool FastAssign = MatrixMajor::isColMatrix<M>();
    };
}
