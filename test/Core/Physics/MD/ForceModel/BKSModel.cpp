/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/BKSModel.h"
#include "Physica/Core/MultiPrecision/AutoDiffGuard.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

using namespace Physica;
using namespace Physica::Core;
using namespace Physica::Utils;

namespace Physica {
    class Test {
        using RandomGenerator = std::mt19937;
        using PlainScalar = Scalar<Double>;
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse, 1>;
        using MDCellType = MDCell<ScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using MassVector = typename MDCellType::MassVector;
        using ForceModel = BKSModel<ScalarType, Ewald<ScalarType>, false>;
        constexpr static size_t MoleculePerCell = 4;
        constexpr static double pair_cutoff = PhyConst<AU>::angstormToBohr(5);
    public:
        static void run() {
            const ScalarType volume = 125;
            const unsigned int cellSize = 3;
            const AutoDiffGuard<ScalarType> guard{};
            RandomGenerator gen{};
            auto cell = makeSystem(cellSize, volume, gen);
            ForceModel forceModel(cell, pair_cutoff, Ewald<ScalarType>{});
            {
                const AutoDiffGuard<ScalarType> guard1{};
                forceModel.potentialEnergy(cell).reverse();
            }
            /* Test press */ {
                const PlainScalar press_diff = -volume.getGrad() / PlainScalar(cellSize * cellSize * cellSize);
                const PlainScalar press = (forceModel.virial(cell).trace() / ScalarType(3)).getValue();
                if (!scalarNear(press_diff, press, 1E-12))
                    exit(EXIT_FAILURE);
            }
        }
    private:
        /**
         * Reference:
         * [1] mp-7000; https://doi.org/10.17188/1272685
         */
        static Vector<ScalarType, 3> randomVector(ScalarType latticeConst, RandomGenerator& gen) {
            constexpr double equalR = PhyConst<AU>::angstormToBohr(1.62844); // Bond length of Si-O, refer to [1]
            std::uniform_real_distribution dist{};
            const ScalarType theta(dist(gen) * M_PI);
            const ScalarType phi(dist(gen) * M_PI * 2);
            Vector<ScalarType, 3> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
            result *= ScalarType(equalR / double(latticeConst.getValue())) * latticeConst; // Pass grad to latticeConst while keep O-H bond length unchanged
            return result;
        }

        static MDCell<ScalarType> makeSystem(unsigned int cellSize, ScalarType cellVolume, RandomGenerator& gen) {
            using CrystalCellType = CrystalCell<ScalarType>;
            using Vector3D = Vector<ScalarType, 3>;
            
            constexpr size_t maxIndexO = MoleculePerCell * 2;
            constexpr size_t maxIndexSi = MoleculePerCell * 3;
            constexpr size_t numAtom = MoleculePerCell * 3;

            const ScalarType latticeConst(cbrt(cellVolume));
            typename CrystalCellType::LatticeMatrix lattice = CrystalCellType::LatticeMatrix::unitMatrix(3);
            lattice *= latticeConst;

            typename CrystalCellType::PositionMatrix pos(numAtom, 3);
            std::uniform_real_distribution dist(-0.1, 0.1);
            for (size_t i = 0; i < MoleculePerCell; ++i) {
                Vector3D temp = Vector3D::random_any(3, dist, gen);
                if (i == 0) {
                    temp[0] += ScalarType(0.25);
                    temp[1] += ScalarType(0.25);
                    temp[2] += ScalarType(0.25);
                }
                else if (i == 1) {
                    temp[0] += ScalarType(0.75);
                    temp[1] += ScalarType(0.75);
                    temp[2] += ScalarType(0.25);
                }
                else if (i == 2) {
                    temp[0] += ScalarType(0.75);
                    temp[1] += ScalarType(0.25);
                    temp[2] += ScalarType(0.75);
                }
                else if (i == 3) {
                    temp[0] += ScalarType(0.25);
                    temp[1] += ScalarType(0.75);
                    temp[2] += ScalarType(0.75);
                }
                if constexpr (ScalarType::isReverseDiff)
                    temp.makeContinuous();
                temp *= latticeConst;
                auto posSi = pos.row(i + maxIndexO);
                auto posO1 = pos.row(2 * i);
                auto posO2 = pos.row(2 * i + 1);
                posSi = temp;
                posO1 = temp + randomVector(latticeConst, gen);
                posO2 = temp + randomVector(latticeConst, gen);
            }

            typename CrystalCellType::AtomicArray atomicNumbers(numAtom);
            for (size_t i = 0; i < maxIndexO; ++i)
                atomicNumbers[i] = 8;
            for (size_t i = maxIndexO; i < maxIndexSi; ++i)
                atomicNumbers[i] = 14;

            CrystalCellType cell({std::move(lattice), std::move(pos), CrystalCellType::Type::Cartesian}, std::move(atomicNumbers));
            cell.toSuperCell(cellSize, cellSize, cellSize);
            cell.normalize();
            return MDCell<ScalarType>(std::move(cell));
        }
    };
}

int main() {
    Physica::Test::run();
    return 0;
}
