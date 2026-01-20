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

#include "Physica/Core/Physics/ManyBody/Model/SquareLattice.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"
#include "HamiltonMatrix.h"

namespace Physica {
    template<Scalar T, int Dim, int NumSite, BoundaryCond BC = BoundaryCond::PBC>
    class J1J2Matrix : public HamiltonMatrix<J1J2Matrix<T, Dim, NumSite, BC>>
                     , public SquareLattice<Dim, BC> {
        using This = J1J2Matrix<T, Dim, NumSite, BC>;
        using Base = HamiltonMatrix<This>;
        using Lattice = SquareLattice<Dim, BC>;
        using Repr = SpinRepr<Dim, NumSite>;
        static_assert(BC == BoundaryCond::PBC, "[Error]: Not implemented");

        using StateType = Repr::StateType;
    public:
        using typename Base::ScalarType;
        constexpr static int SiteDOF = StateType::SiteDOF;
    protected:
        using typename Base::Tr;
        using typename Base::Trv;
    private:
        Tr couplingJ1;
        Tr couplingJ2;
        Repr repr;
    public:
        J1J2Matrix() = default;
        J1J2Matrix(Tr couplingJ1, Tr couplingJ2, Lattice lattice);
        J1J2Matrix(const This&) = default;
        J1J2Matrix(This&&) noexcept = default;
        ~J1J2Matrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Tr getCouplingJ1() const noexcept { return couplingJ1; }
        [[nodiscard]] Tr getCouplingJ2() const noexcept { return couplingJ2; }
        [[nodiscard]] const Lattice& getModel() const noexcept { return *this; }
        [[nodiscard]] const Repr& getRepr() const noexcept { return repr; }
    };

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    J1J2Matrix<T, Dim, NumSite, BC>::J1J2Matrix(Tr couplingJ1, Tr couplingJ2, Lattice lattice)
            : Lattice(std::move(lattice)), couplingJ1(couplingJ1), couplingJ2(couplingJ2) {
        assert(Lattice::getNumSuperCellSite() == NumSite && "[Error]: Number is not consistent");
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    auto J1J2Matrix<T, Dim, NumSite, BC>::operator*(this auto&& self, Vector auto&& v) noexcept {
        using Self = decltype(self);
        using V = decltype(v);
        return GEMV<Self, V>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    T J1J2Matrix<T, Dim, NumSite, BC>::calc(size_t row, size_t col) const {
        using IntType = StateType::IntType;
        const auto psiR = repr[row];
        if (row == col) {
            int numAntiSpin = 0;
            int numNextAntiSpin = 0;
            for (int i = 0; i < NumSite; ++i) {
                Lattice::forNeighSites([&, psiR](int from, int to) noexcept {
                    numAntiSpin += psiR[from] == psiR[to] ? 1 : -1;
                }, i);

                Lattice::forNNeighSites([&, psiR](int from, int to) noexcept {
                    numNextAntiSpin += psiR[from] == psiR[to] ? 1 : -1;
                }, i);
            }
            return (getCouplingJ1() * T(numAntiSpin) + getCouplingJ2() * T(numNextAntiSpin)) * Trv(0.25);
        }

        const auto psiC = repr[col];
        const auto diffMask = psiR.getOccupyBits() ^ psiC.getOccupyBits();
        bool isParallel = std::popcount(psiR.getOccupyBits() & diffMask) != 1;
        int numSpinDiff = std::popcount(diffMask);
        if ((numSpinDiff != 2) || isParallel)
            return Trv(0);

        bool match = false;
        auto finder = [&, diffMask](int from, int to) noexcept {
            IntType mask = (IntType(1) << from) | (IntType(1) << to);
            match |= mask == diffMask;
        };

        for (int i = 0; i < NumSite; ++i) {
            Lattice::forNeighSites(finder, i);
            if (match)
                return Trv(0.5) * getCouplingJ1();

            Lattice::forNNeighSites(finder, i);
            if (match)
                return Trv(0.5) * getCouplingJ2();
        }
        return Trv(0);
    }

    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    void J1J2Matrix<T, Dim, NumSite, BC>::swap(This& __restrict obj) noexcept {
        Lattice::swap(obj);
        couplingJ1.swap(obj.couplingJ1);
        couplingJ2.swap(obj.couplingJ2);
        repr.swap(obj.repr);
    }
}

namespace Physica {
    template<Scalar T, int Dim, int NumSite, BoundaryCond BC>
    class Traits<J1J2Matrix<T, Dim, NumSite, BC>> : public Traits<HamiltonMatrix<J1J2Matrix<T, Dim, NumSite, BC>>> {
    public:
        using ScalarType = T;
        using ReprType = SpinRepr<Dim, NumSite>;
        constexpr static BoundaryCond Boundary = BC;
    };
}

#include "J1J2GEMV.h"
