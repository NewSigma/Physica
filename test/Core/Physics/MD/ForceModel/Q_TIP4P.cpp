/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"

using namespace Physica;
using RandomSource = Random<MT19937, 3251508104675666787>;

// FIXME: Re-enable this test
constexpr bool Disable = sizeof(int) == 0;

namespace Physica {
    class Test {
        using T = Diff<float64, DiffMode::Reverse, 1>;
        using Tv = T::ValueType;
        using MDCellType = MDCell<T>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using MassVector = MDCellType::MassVector;
        using ForceModel = Q_TIP4P<T, Ewald<T>>;
        constexpr static size_t MoleculePerCell = 4;
        constexpr static double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
        constexpr static double massMoleculeInSI = (PhyConst<SI>::atomMass(1) * 2) + PhyConst<SI>::atomMass(8);
    public:
        static void run() {
            if constexpr (Disable) {
                const Tv volume = ((MoleculePerCell * massMoleculeInSI * 1000 / 0.997) * 1E-6) / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius);
                const unsigned int cellSize = 3;
                auto cell = makeSystem(cellSize, volume);
                ForceModel::sortPosition(cell);
                ForceModel forceModel(cell, pair_cutoff, Ewald<T>{});
                forceModel.potentialV(cell).reverse();
                /* Test press */ {
                    const Tv press_diff = -volume.grad() / Tv(cellSize * cellSize * cellSize);
                    const Tv press = (forceModel.virial(cell).trace() / Tv(3)).value();
                    expect(scalarNear(press_diff, press, 1E-12));
                }
            }
        }
    private:
        static Vector3D<T> randomVector(T latticeConst) {
            const Tv theta(Tv::random_uniform<RandomSource>() * Tv(M_PI));
            const Tv phi(Tv::random_uniform<RandomSource>() * Tv(M_PI * 2));
            Vector3D<T> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
            result *= T(ForceModel::equalR / double(latticeConst.value())) * latticeConst; // Pass grad to latticeConst while keep O-H bond length unchanged
            return result;
        }

        static MDCell<T> makeSystem(unsigned int cellSize, Tv cellVolume) {
            if constexpr (Disable) {
                using CrystalCellType = CrystalCell<T>;
                constexpr size_t maxIndexH = MoleculePerCell * 2;
                constexpr size_t maxIndexO = MoleculePerCell * 3;
                constexpr size_t numAtom = MoleculePerCell * 3;

                const Tv latticeConst(cbrt(cellVolume));
                CrystalCellType::LatticeMatrix lattice = CrystalCellType::LatticeMatrix::identity(3);
                lattice *= latticeConst;

                CrystalCellType::PositionMatrix pos(numAtom, 3);
                std::uniform_real_distribution dist(-0.1, 0.1);
                for (size_t i = 0; i < MoleculePerCell; ++i) {
                    auto temp = Vector3D<Tv>::random_any<RandomSource>(3, dist);
                    if (i == 0) {
                        temp[0] += Tv(0.25);
                        temp[1] += Tv(0.25);
                        temp[2] += Tv(0.25);
                    }
                    else if (i == 1) {
                        temp[0] += Tv(0.75);
                        temp[1] += Tv(0.75);
                        temp[2] += Tv(0.25);
                    }
                    else if (i == 2) {
                        temp[0] += Tv(0.75);
                        temp[1] += Tv(0.25);
                        temp[2] += Tv(0.75);
                    }
                    else if (i == 3) {
                        temp[0] += Tv(0.25);
                        temp[1] += Tv(0.75);
                        temp[2] += Tv(0.75);
                    }

                    temp *= latticeConst;
                    auto posO = pos.row(i + maxIndexH);
                    auto posH1 = pos.row(2 * i);
                    auto posH2 = pos.row(2 * i + 1);
                    posO = temp;
                    posH1 = temp + randomVector(latticeConst);
                    posH2 = temp + randomVector(latticeConst);
                }

                CrystalCellType::AtomicArray atomicNumbers(numAtom);
                for (size_t i = 0; i < maxIndexH; ++i)
                    atomicNumbers[i] = 1;
                for (size_t i = maxIndexH; i < maxIndexO; ++i)
                    atomicNumbers[i] = 8;

                CrystalCellType cell({std::move(lattice), std::move(pos), CrystalCellType::Type::Cartesian}, std::move(atomicNumbers));
                cell.toSuperCell(cellSize, cellSize, cellSize);
                return MDCell<T>(std::move(cell));
            }
            return {};
        }
    };
}

int main() {
    Physica::Test::run();
    return 0;
}
