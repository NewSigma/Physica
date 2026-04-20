/*
 * Copyright 2025-2026 Weibo He.
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

#include "AndersonMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<AndersonMatrix, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        constexpr static int NumSite = Traits<std::remove_cvref_t<M>>::NumSite;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M mat_, V vec_);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const noexcept;

        [[nodiscard]] T calc(size_t) const { noImpl(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M, Vector V> requires(instanceof_tx<AndersonMatrix, M>)
    GEMV<M, V>::GEMV(M mat_, V vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {}

    template<Matrix M, Vector V> requires(instanceof_tx<AndersonMatrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const noexcept {
        target.assert_assign(*this);
        target.zeros();

        const auto& repr = mat.getRepr();
        for (size_t i = 0; i < target.getLength(); ++i) {
            const auto psiC = repr[i];
            for (int site = 0; site < NumSite - 1; ++site) {
                int site1 = site + 1;
                const bool upOccupy1 = psiC.isUpOccupy(site);
                const bool downOccupy1 = psiC.isDownOccupy(site);
                const bool upOccupy2 = psiC.isUpOccupy(site1);
                const bool downOccupy2 = psiC.isDownOccupy(site1);
                const T hop = mat.getHoppingT(site) * vec.calc(i);
                if (upOccupy1 != upOccupy2) {
                    const auto psiR = upOccupy1 ? psiC.hopUp(site, site1) : psiC.hopUp(site1, site);
                    const size_t index = repr[psiR];
                    const bool sign = upOccupy1 == (psiC.hopUpSign(site, site1) == 1);
                    target[index] += (sign ? hop : -hop);
                }

                if (downOccupy1 != downOccupy2) {
                    const auto psiR = downOccupy1 ? psiC.hopDown(site, site1) : psiC.hopDown(site1, site);
                    const size_t index = repr[psiR];
                    const bool sign = downOccupy1 == (psiC.hopDownSign(site, site1) == 1);
                    target[index] += (sign ? hop : -hop);
                }
            }
            target[i] += mat.calc(i, i) * vec.calc(i);
        }
    }

    template<Matrix M, Vector V> requires(instanceof_tx<AndersonMatrix, M>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<AndersonMatrix, M>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<AndersonMatrix, M>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref<M>::type;
        using V1 = std::remove_cvref<V>::type;
        using T1 = M1::ScalarType;
        using T2 = V1::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        constexpr static bool FastAssign = true;
    };
}
