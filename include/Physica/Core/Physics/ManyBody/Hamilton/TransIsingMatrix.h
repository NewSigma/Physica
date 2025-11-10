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

#include "Physica/Core/Physics/ManyBody/Model/SquareLattice.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"
#include "HamiltonMatrix.h"

namespace Physica {
    template<Scalar, int, int, BoundaryCond, Vector>
    class TransIsingVecProd;

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC = BoundaryCond::PBC>
    class TransIsingMatrix : public HamiltonMatrix<TransIsingMatrix<T, Dim, NumSite, BC>>
                           , public SquareLattice<Dim, BC> {
        using This = TransIsingMatrix<T, Dim, NumSite, BC>;
        using Base = HamiltonMatrix<This>;
        using Lattice = SquareLattice<Dim, BC>;
        using Repr = SpinRepr<Dim, NumSite>;

        using StateType = Repr::StateType;
    public:
        using typename Base::ScalarType;
        constexpr static int SiteDOF = StateType::SiteDOF;
    private:
        T couplingJ;
        T transG;
        Repr repr;
    public:
        TransIsingMatrix() = default;
        TransIsingMatrix(T couplingJ_, T transG_, Lattice lattice);
        TransIsingMatrix(const This&) = default;
        TransIsingMatrix(This&&) noexcept = default;
        ~TransIsingMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        [[nodiscard]] auto operator*(const Vector auto& v) const noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getCouplingJ() const noexcept { return couplingJ; }
        [[nodiscard]] T getTransG() const noexcept { return transG; }
        [[nodiscard]] const Lattice& getModel() const noexcept { return *this; }
        [[nodiscard]] const Repr& getRepr() const noexcept { return repr; }
    };

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    TransIsingMatrix<T, Dim, NumSite, BC>::TransIsingMatrix(T couplingJ_, T transG_, Lattice lattice)
            : Lattice(std::move(lattice)), couplingJ(couplingJ_), transG(transG_) {
        assert(Lattice::getNumSuperCellSite() == repr.getNumSpin() && "[Error]: Number is not consistent");
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    auto TransIsingMatrix<T, Dim, NumSite, BC>::operator*(const Vector auto& v) const noexcept {
        using V = std::remove_cvref_t<decltype(v)>;
        return TransIsingVecProd<T, Dim, NumSite, BC, V>(*this, v);
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    T TransIsingMatrix<T, Dim, NumSite, BC>::calc(size_t row, size_t col) const {
        const auto psi1 = repr[row];
        if (row == col) {
            using enum BoundaryCond;
            const auto bits = psi1.getOccupyBits();
            if constexpr (BC == OBC) {
                constexpr auto Mask = StateType::makeFullMask() >> 1;
                const int numAntiSpin = std::popcount((bits >> 1) ^ (bits & Mask));
                assert(numAntiSpin < NumSite);
                return -getCouplingJ() * T(NumSite - 1 - numAntiSpin * 2);
            }
            else {
                static_assert(BC == PBC, "[Error]: Unsupported boundary condition");
                constexpr auto Mask = StateType::makeHighMask();
                const int numAntiSpin = std::popcount(((bits >> 1) | ((bits & 1) ? Mask : 0)) ^ bits);
                return -getCouplingJ() * T(NumSite - numAntiSpin * 2);
            }
        }

        const auto psi2 = repr[col];
        bool OneFlip = (psi1 ^ psi2).getNumParticle() == 1;
        if (OneFlip)
            return -getTransG();
        return 0;
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    void TransIsingMatrix<T, Dim, NumSite, BC>::swap(This& __restrict obj) noexcept {
        Lattice::swap(obj);
        couplingJ.swap(obj.couplingJ);
        transG.swap(obj.transG);
        repr.swap(obj.repr);
    }
}

namespace Physica {
    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    class Traits<TransIsingMatrix<T, Dim, NumSite, BC>> : public Traits<HamiltonMatrix<TransIsingMatrix<T, Dim, NumSite, BC>>> {
    public:
        using ScalarType = T;
        using ReprType = SpinRepr<Dim, NumSite>;
        constexpr static BoundaryCond Boundary = BC;
    };
}

#include "TransIsingVecProd.h" // IWYU pragma: export
