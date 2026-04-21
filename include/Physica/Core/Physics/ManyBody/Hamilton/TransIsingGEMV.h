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

#include "Physica/Core/Physics/ManyBody/ReprSpace/PauliMatrix.h"
#include "TransIsingMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<TransIsingMatrix, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        using M1 = std::remove_cvref<M>::type;
        constexpr static unsigned int NumSite = M1::NumSite;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M&& mat_, V&& vec_);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] T calc(size_t) const { noImpl(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    private:
        /* Getters */
        [[nodiscard]] const auto& getRepr() const noexcept { return mat.getRepr(); }
        [[nodiscard]] T getCouplingJ() const noexcept { return mat.getCouplingJ(); }
        [[nodiscard]] T getTransG() const noexcept { return mat.getTransG(); }
    };

    template<Matrix M, Vector V> requires(instanceof_tx<TransIsingMatrix, M>)
    GEMV<M, V>::GEMV(M&& mat_, V&& vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_tx<TransIsingMatrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        target.assert_assign(*this);
        target.zeros();

        const auto& repr = getRepr();
        for (size_t i = 0; i < getLength(); ++i) {
            using enum PauliIndex;
            auto psi0 = repr[i];
            const T g = getTransG() * vec.calc(i);
            for (unsigned int site = 0; site < NumSite; ++site) {
                auto psi = psi0;
                target[repr[psi]] -= PauliMatrix<T, X>(site).apply(psi) * g;
            }

            const T couplingJ = getCouplingJ() * vec.calc(i);
            constexpr BoundaryCond BC = Traits<M>::Boundary;
            if constexpr (BC == BoundaryCond::OBC) {
                for (unsigned int site = 0; site < NumSite - 1; ++site) {
                    auto coeff = PauliMatrix<T, Z>(site).apply(psi0);
                    coeff *= PauliMatrix<T, Z>(site + 1).apply(psi0);
                    target[i] -= couplingJ * coeff;
                }
            }
            else {
                static_assert(BC == BoundaryCond::PBC, "[Error]: Unsupported BoundaryCond");
                for (unsigned int site = 0; site < NumSite; ++site) {
                    auto coeff = PauliMatrix<T, Z>(site).apply(psi0);
                    coeff *= PauliMatrix<T, Z>((site + 1) % NumSite).apply(psi0);
                    target[i] -= couplingJ * coeff;
                }
            }
        }
    }

    template<Matrix M, Vector V> requires(instanceof_tx<TransIsingMatrix, M>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<TransIsingMatrix, M>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<TransIsingMatrix, M>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref<M>::type;
        using V1 = std::remove_cvref<V>::type;
        using T1 = M1::ScalarType;
        using T2 = V1::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
    };
}
