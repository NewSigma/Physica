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
#include "Physica/Core/Physics/MD/ForceModel/BKSModel.h"

using namespace Physica;
using namespace Physica;

namespace Physica {
    class Test {
        using RandomType = Random<MT19937>;
        using dfloat = Diff<float64, DiffMode::Reverse, 1>;
        using MDCellType = MDCell<dfloat>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using MassVector = MDCellType::MassVector;
        using ForceModel = BKSModel<dfloat, Ewald<dfloat>, false>;
        constexpr static size_t MoleculePerCell = 4;
        constexpr static double pair_cutoff = PhyConst<AU>::angstormToBohr(5);
    public:
        static void run() {
            dfloat volume = 125;
            const unsigned int cellSize = 3;
            auto cell = makeSystem(volume, cellSize);
            ForceModel forceModel(cell, pair_cutoff, Ewald<dfloat>{});
            forceModel.potentialV(cell).reverse();
            /* Test press */ {
                const float64 press_diff = -volume.grad() / float64(cellSize * cellSize * cellSize);
                const float64 press = (forceModel.virial(cell).trace() / float64(3)).value();
                if (!scalarNear(press_diff, press, 1E-12))
                    exit(EXIT_FAILURE);
            }
        }
    private:
        /**
         * Reference:
         * [1] mp-7000; https://doi.org/10.17188/1272685
         */
        static Vector3D<float64> randomVector(float64 latticeConst) {
            constexpr double equalR = PhyConst<AU>::angstormToBohr(1.62844); // Bond length of Si-O, refer to [1]
            std::uniform_real_distribution dist{};
            const auto theta = float64::random_uniform<RandomType>() * float64(M_PI);
            const auto phi = float64::random_uniform<RandomType>() * float64(M_PI * 2);
            Vector3D<float64> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
            result *= float64(equalR / double(latticeConst.value())) * latticeConst; // Pass grad to latticeConst while keep O-H bond length unchanged
            return result;
        }

        static MDCell<dfloat> makeSystem(dfloat& cellVolume, unsigned int cellSize) {
            using CrystalCellType = CrystalCell<float64>;
            using LatticeMatrix = CrystalCellType::LatticeMatrix;
            constexpr size_t maxIndexO = MoleculePerCell * 2;
            constexpr size_t maxIndexSi = MoleculePerCell * 3;
            constexpr size_t numAtom = MoleculePerCell * 3;

            const dfloat latticeConst(cbrt(cellVolume));
            CoDiff<LatticeMatrix> lattice = LatticeMatrix::unitMatrix(3) * latticeConst;

            CrystalCellType::PositionMatrix pos(numAtom, 3);
            std::uniform_real_distribution dist(-0.1, 0.1);
            for (size_t i = 0; i < MoleculePerCell; ++i) {
                auto temp = Vector3D<float64>::template random_any<R, decltype(dist)>(3, dist);
                if (i == 0) {
                    temp[0] += float64(0.25);
                    temp[1] += float64(0.25);
                    temp[2] += float64(0.25);
                }
                else if (i == 1) {
                    temp[0] += float64(0.75);
                    temp[1] += float64(0.75);
                    temp[2] += float64(0.25);
                }
                else if (i == 2) {
                    temp[0] += float64(0.75);
                    temp[1] += float64(0.25);
                    temp[2] += float64(0.75);
                }
                else if (i == 3) {
                    temp[0] += float64(0.25);
                    temp[1] += float64(0.75);
                    temp[2] += float64(0.75);
                }

                temp *= latticeConst;
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
    };
}

int main() {
    Physica::Test::run();
    return 0;
}
