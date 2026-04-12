/*
 * Copyright 2025 Weibo He.
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
#include "J1J2Matrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<J1J2Matrix, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        using M1 = std::remove_cvref<M>::type;
        constexpr static unsigned int NumSite = M1::NumSite;
    protected:
        using typename Base::T;
        using typename Base::Trv;
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
    private:
        /* Getters */
        [[nodiscard]] const auto& getRepr() const noexcept { return mat.getRepr(); }
        [[nodiscard]] auto getCouplingJ1() const noexcept { return mat.getCouplingJ1(); }
        [[nodiscard]] auto getCouplingJ2() const noexcept { return mat.getCouplingJ2(); }
    };

    template<Matrix M, Vector V> requires(instanceof_tx<J1J2Matrix, M>)
    GEMV<M, V>::GEMV(M&& mat_, V&& vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_tx<J1J2Matrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        target.assert_assign(*this);
        target.zeros();

        const auto& repr = getRepr();
        for (size_t i = 0; i < getLength(); ++i) {
            using enum PauliIndex;
            auto psiR = repr[i];
            /* Z component */ {
                int numAntiSpin = 0;
                int numNextAntiSpin = 0;
                for (int site = 0; site < NumSite; ++site) {
                    mat.forNeighSites([&, psiR](int from, int to) noexcept {
                        numAntiSpin += psiR[from] == psiR[to] ? 1 : -1;
                    }, site);

                    mat.forNNeighSites([&, psiR](int from, int to) noexcept {
                        numNextAntiSpin += psiR[from] == psiR[to] ? 1 : -1;
                    }, site);
                }
                target[i] += Trv(0.25) * vec.calc(i) * (getCouplingJ1() * T(numAntiSpin) + getCouplingJ2() * T(numNextAntiSpin));
            }
            /* XY component */ {
                for (int site = 0; site < NumSite; ++site) {
                    mat.forNeighSites([&, psiR](int from, int to) noexcept {
                        if (psiR[from] != psiR[to]) {
                            auto psiC = psiR.flip(from).flip(to);
                            target[repr[psiC]] += vec.calc(i) * getCouplingJ1() * T(0.5);
                        }
                    }, site);

                    mat.forNNeighSites([&, psiR](int from, int to) noexcept {
                        if (psiR[from] != psiR[to]) {
                            auto psiC = psiR.flip(from).flip(to);
                            target[repr[psiC]] += vec.calc(i) * getCouplingJ2() * T(0.5);
                        }
                    }, site);
                }
            }
        }
    }

    template<Matrix M, Vector V> requires(instanceof_tx<J1J2Matrix, M>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<J1J2Matrix, M>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<J1J2Matrix, M>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref<M>::type;
        using V1 = std::remove_cvref<V>::type;
        using T1 = M1::ScalarType;
        using T2 = V1::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        constexpr static size_t SizeAtCompile = M1::RowAtCompile;
        constexpr static bool FastAssign = true;
    };
}
