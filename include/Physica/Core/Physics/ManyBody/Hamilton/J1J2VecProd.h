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
    template<Scalar U, int Dim, int NumSite, BoundaryCond BC, Vector V>
    class J1J2VecProd : public RValueVector<J1J2VecProd<U, Dim, NumSite, BC, V>> {
        using MatrixType = J1J2Matrix<U, Dim, NumSite, BC>;
        using This = J1J2VecProd<U, Dim, NumSite, BC, V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Trv;
    private:
        const MatrixType& mat;
        const V& vec;
    public:
        J1J2VecProd(const MatrixType& mat_, const V& vec_);
        J1J2VecProd(const This&) = default;
        J1J2VecProd(This&&) noexcept = default;
        ~J1J2VecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] T calc(size_t) const { noImpl(__func__); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const V& getRHS() const noexcept { return vec; }
    private:
        /* Getters */
        [[nodiscard]] const auto& getRepr() const noexcept { return mat.getRepr(); }
        [[nodiscard]] auto getCouplingJ1() const noexcept { return mat.getCouplingJ1(); }
        [[nodiscard]] auto getCouplingJ2() const noexcept { return mat.getCouplingJ2(); }
    };

    template<Scalar U, int Dim, int NumSite, BoundaryCond BC, Vector V>
    J1J2VecProd<U, Dim, NumSite, BC, V>::J1J2VecProd(const MatrixType& mat_, const V& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Scalar U, int Dim, int NumSite, BoundaryCond BC, Vector V>
    template<ExecutePolicy P>
    void J1J2VecProd<U, Dim, NumSite, BC, V>::assign(Vector auto& target) const {
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
}

namespace Physica {
    template<Scalar U, int Dim, int NumSite, BoundaryCond BC, Vector V>
    class Traits<J1J2VecProd<U, Dim, NumSite, BC, V>> {
        using MatrixType = J1J2Matrix<U, Dim, NumSite, BC>;
        using U1 = V::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<U, U1>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
