/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Test.h"

using namespace Physica;

namespace {
    template<Scalar T>
    void VASPTest() {
        constexpr static bool isFloat32 = T::Prec == Float32;
        constexpr double prec = isFloat32 ? 2E-5 : 1E-5;

        const double lengthInBohr = PhyConst<AU>::angstromToBohr(3);
        CrystalCell<T> cell({{lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {0.5, 0.5, 0.5}, CrystalCell<T>::Type::Direct}, {14});
        Ewald<T> ewald(cell.getLattice(), {4});
        const auto energy = ewald.potentialV(cell.getPos());
        expect(scalarNear(energy, T(PhyConst<AU>::eVToHartree(-108.95061336198556)), prec));
    }
    /**
    * Reference:
    * [1] pyewald; https://github.com/lukeolson/pyewald
    */
    template<Scalar T>
    void madelungTest() {
        constexpr static bool isFloat32 = T::Prec == Float32;
        using EwaldType = Ewald<T>;
        {
            const double lengthInBohr = PhyConst<AU>::angstromToBohr(5.6903014761756712);
            CrystalCell<T> NaCl({{lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {
                                0.0, 0.0, 0.0,
                                0.0, 0.5, 0.5,
                                0.5, 0.0, 0.5,
                                0.5, 0.5, 0.0,
                                0.5, 0.5, 0.5,
                                0.5, 0.0, 0.0,
                                0.0, 0.5, 0.0,
                                0.0, 0.0, 0.5
                            }, CrystalCell<T>::Type::Direct}, {1, 1, 1, 1, 1, 1, 1, 1});
            NaCl.toCartesian();
            EwaldType ewald(NaCl.getLattice(), {1, 1, 1, 1, -1, -1, -1, -1});
            const auto energy = ewald.potentialV(NaCl.getPos());
            const auto madelung = -(energy / 4) * (lengthInBohr / 2); // We have 4x  cell so energy is divided by 4
            constexpr double prec = isFloat32 ? 1E-6 : 1E-7;
            expect(scalarNear(madelung, T(1.7475645946331822), prec));
        }
        {
            const double lengthInBohr = 1;
            CrystalCell<T> CsCl({{lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {
                                0.0, 0.0, 0.0,
                                0.5, 0.5, 0.5,
                            }, CrystalCell<T>::Type::Cartesian}, {1, 1});
            EwaldType ewald(CsCl.getLattice(), {1, -1});
            const auto energy = ewald.potentialV(CsCl.getPos());
            const auto madelung = -energy * (lengthInBohr * 0.5 * std::sqrt(3.0));
            constexpr double prec = isFloat32 ? 1E-5 : 1E-9;
            expect(scalarNear(madelung, T(1.76267477307099), prec));
        }
        {
            const double lengthInBohr = 0.5;
            CrystalCell<T> ZnS({{0, lengthInBohr, lengthInBohr, lengthInBohr, 0, lengthInBohr, lengthInBohr, lengthInBohr, 0}, {
                                0.0, 0.0, 0.0,
                                0.25, 0.25, 0.25,
                            }, CrystalCell<T>::Type::Direct}, {1, 1});
            ZnS.toCartesian();
            EwaldType ewald(ZnS.getLattice(), {1, -1});
            const auto energy = ewald.potentialV(ZnS.getPos());
            const auto madelung = -energy * (lengthInBohr * 0.5 * std::sqrt(3.0));
            constexpr double prec = isFloat32 ? 1E-5 : 1E-9;
            expect(scalarNear(madelung, T(1.63805505338879), prec));
        }
    }
}

int main() {
    VASPTest<float64>();
    VASPTest<float32>();
    madelungTest<float64>();
    madelungTest<float32>();
    return 0;
}
