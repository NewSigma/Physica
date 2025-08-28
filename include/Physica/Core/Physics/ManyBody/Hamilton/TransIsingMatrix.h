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

#include "Physica/Core/Physics/ManyBody/ReprSpace/ReprBasis.h"
#include "Physica/Core/Physics/ManyBody/Model/SquareLattice.h"
#include "HamiltonMatrix.h"

namespace Physica {
    template<Scalar, Representation, BoundaryCond, Vector>
    class TransIsingVecProd;

    template<Scalar T, Representation Repr, BoundaryCond BC = BoundaryCond::PBC>
    class TransIsingMatrix : public HamiltonMatrix<TransIsingMatrix<T, Repr, BC>>
                           , public SquareLattice<Repr::Dim, BC> {
        using This = TransIsingMatrix<T, Repr, BC>;
        using Base = HamiltonMatrix<This>;
        using Lattice = SquareLattice<Repr::Dim, BC>;

        using StateType = Repr::StateType;
    public:
        using typename Base::ScalarType;
        constexpr static int Dim = Repr::Dim;
        constexpr static int NumSite = StateType::NumSite;
        constexpr static int SiteDOF = StateType::SiteDOF;
    private:
        T couplingJ;
        T transG;
        Repr repr;
    public:
        TransIsingMatrix() = default;
        TransIsingMatrix(T couplingJ_, T transG_, Lattice lattice, Repr repr_);
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

    template<Scalar T, Representation Repr, BoundaryCond BC>
    TransIsingMatrix<T, Repr, BC>::TransIsingMatrix(T couplingJ_, T transG_, Lattice lattice, Repr repr_)
            : Lattice(std::move(lattice)), couplingJ(couplingJ_), transG(transG_), repr(std::move(repr_)) {}

    template<Scalar T, Representation Repr, BoundaryCond BC>
    auto TransIsingMatrix<T, Repr, BC>::operator*(const Vector auto& v) const noexcept {
        using V = std::remove_cvref_t<decltype(v)>;
        return TransIsingVecProd<T, Repr, BC, V>(*this, v);
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    T TransIsingMatrix<T, Repr, BC>::calc(size_t row, size_t col) const {
        const auto psi1 = repr[row];
        const auto mask = StateType::makeFullMask() >> 1;
        if (row == col) {
            const auto bits = psi1.getOccupyBits();
            const int numAntiSpin = std::popcount((bits >> 1) ^ (bits & mask));
            assert(numAntiSpin < NumSite);
            return -getCouplingJ() * T((int(NumSite - 1) - numAntiSpin * 2));
        }

        const auto psi2 = repr[col];
        if ((psi1 ^ psi2).getNumParticle() == 1)
            return -getTransG();
        return 0;
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    void TransIsingMatrix<T, Repr, BC>::swap(This& __restrict obj) noexcept {
        Lattice::swap(obj);
        couplingJ.swap(obj.couplingJ);
        transG.swap(obj.transG);
        repr.swap(obj.repr);
    }
}

namespace Physica {
    template<Scalar T, Representation Repr, BoundaryCond BC>
    class Traits<TransIsingMatrix<T, Repr, BC>> : public Traits<HamiltonMatrix<TransIsingMatrix<T, Repr, BC>>> {
    public:
        using ScalarType = T;
        using ReprType = Repr;
        constexpr static BoundaryCond Boundary = BC;
    };
}

#include "TransIsingVecProd.h" // IWYU pragma: export
