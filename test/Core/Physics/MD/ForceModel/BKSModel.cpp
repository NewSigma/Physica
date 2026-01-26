/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/BKSModel.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

// FIXME: Re-enable this test
constexpr bool Disable = sizeof(int) == 0;

namespace Physica {
    class Test {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Reverse, 1>;
        using MDCellType = MDCell<dfloat>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using MassVector = MDCellType::MassVector;
        using ForceModel = BKSModel<dfloat, Ewald<dfloat>, false>;
        constexpr static size_t MoleculePerCell = 4;
        constexpr static double pair_cutoff = PhyConst<AU>::angstormToBohr(5);
    public:
        static void run() {
            if constexpr (Disable) {
                dfloat volume = 125;
                const unsigned int cellSize = 3;
                auto cell = makeSystem(volume, cellSize);
                ForceModel forceModel(cell, pair_cutoff, Ewald<dfloat>{});
                forceModel.potentialV(cell).reverse();
                /* Test press */ {
                    const T press_diff = -volume.grad() / T(cellSize * cellSize * cellSize);
                    const T press = (forceModel.virial(cell).trace() / T(3)).value();
                    expect(scalarNear(press_diff, press, 1E-12));
                }
            }
        }
    private:
        /**
         * Reference:
         * [1] mp-7000; https://doi.org/10.17188/1272685
         */
        static Vector3D<T> randomVector(T latticeConst) {
            constexpr double equalR = PhyConst<AU>::angstormToBohr(1.62844); // Bond length of Si-O, refer to [1]
            const auto theta = T::random_uniform<RandomSource>() * T(M_PI);
            const auto phi = T::random_uniform<RandomSource>() * T(M_PI * 2);
            Vector3D<T> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
            result *= T(equalR / double(latticeConst.value())) * latticeConst; // Pass grad to latticeConst while keep O-H bond length unchanged
            return result;
        }

        static MDCell<dfloat> makeSystem(dfloat& cellVolume, unsigned int cellSize) {
            if constexpr (Disable) {
                using CrystalCellType = CrystalCell<T>;
                using LatticeMatrix = CrystalCellType::LatticeMatrix;
                constexpr size_t maxIndexO = MoleculePerCell * 2;
                constexpr size_t maxIndexSi = MoleculePerCell * 3;
                constexpr size_t numAtom = MoleculePerCell * 3;

                const auto latticeConst = cbrt(cellVolume);
                CoDiff<LatticeMatrix> lattice = LatticeMatrix::identity(3) * latticeConst;

                CrystalCellType::PositionMatrix pos(numAtom, 3);
                std::uniform_real_distribution dist(-0.1, 0.1);
                for (size_t i = 0; i < MoleculePerCell; ++i) {
                    auto temp = Vector3D<T>::template random_any<RandomSource>(3, dist);
                    if (i == 0) {
                        temp[0] += T(0.25);
                        temp[1] += T(0.25);
                        temp[2] += T(0.25);
                    }
                    else if (i == 1) {
                        temp[0] += T(0.75);
                        temp[1] += T(0.75);
                        temp[2] += T(0.25);
                    }
                    else if (i == 2) {
                        temp[0] += T(0.75);
                        temp[1] += T(0.25);
                        temp[2] += T(0.75);
                    }
                    else if (i == 3) {
                        temp[0] += T(0.25);
                        temp[1] += T(0.75);
                        temp[2] += T(0.75);
                    }

                    temp *= latticeConst.value();
                    auto posSi = pos.row(i + maxIndexO);
                    auto posO1 = pos.row(2 * i);
                    auto posO2 = pos.row(2 * i + 1);
                    posSi = temp;
                    posO1 = temp + randomVector(latticeConst.value());
                    posO2 = temp + randomVector(latticeConst.value());
                }

                CrystalCellType::AtomicArray atomicNumbers(numAtom);
                for (size_t i = 0; i < maxIndexO; ++i)
                    atomicNumbers[i] = 8;
                for (size_t i = maxIndexO; i < maxIndexSi; ++i)
                    atomicNumbers[i] = 14;

                CrystalCellType cell({std::move(lattice), std::move(pos), CrystalCellType::Type::Cartesian}, std::move(atomicNumbers));
                cell.toSuperCell(cellSize, cellSize, cellSize);
                cell.normalize();
                return MDCell<dfloat>(std::move(cell));
            }
            return {};
        }
    };
}

int main() {
    Physica::Test::run();
    return 0;
}
