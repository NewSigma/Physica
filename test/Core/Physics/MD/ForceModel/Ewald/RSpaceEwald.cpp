/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/RSpaceEwald.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/MultiPrecision/AutoDiffGuard.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

using namespace Physica;
using namespace Physica::Core;

namespace Physica {
    class Test {
        using ValueType = float64;
        using ScalarType = Diff<ValueType, DiffMode::Reverse, 1>;
        using MDCellType = MDCell<ScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using MassVector = typename MDCellType::MassVector;
    public:
        static void run() {
            const ScalarType volume = 125;
            const size_t cellSize = 3;
            const AutoDiffGuard<ScalarType> guard{};
            const auto cell = makeSystem(volume, cellSize);
            const auto& pos = cell.getPos();
            VectorND<ScalarType> charges(cell.getNumParticle(), 1.0);
            auto tail = charges.tail(cell.getNumParticle() / 2);
            tail = ScalarType(-1);
            RSpaceEwald<ScalarType, false> ewald(cell.getLattice(), std::move(charges));
            {
                const AutoDiffGuard<ScalarType> guard1{};
                ewald.potentialV(pos).reverse();
            }
            /* Test press */ {
                const ValueType press_diff = -volume.getGrad() / ValueType(cellSize * cellSize * cellSize);
                const ValueType press = (ewald.virial(pos).trace() / ScalarType(3)).getValue();
                if (!scalarNear(press_diff, press, 1E-13))
                    exit(EXIT_FAILURE);
            }
        }
    private:
        static MDCellType makeSystem(ScalarType volume, size_t cellSize) {
            constexpr size_t numMolecularUnitCell = 2;
            LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
            const ScalarType latticeConst = cbrt(volume);
            lattice *= latticeConst;

            PositionMatrix pos(numMolecularUnitCell, 3, ValueType(0));
            auto row = pos.row(1);
            row = ScalarType(0.5);
            pos *= latticeConst;

            MassVector massVec(numMolecularUnitCell, 1);
            MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));
            cell.toSuperCell<ExtendCellOption::AtomMajor>({cellSize, cellSize, cellSize});
            return cell;
        }
    };
}

int main() {
    Physica::Test::run();
    return 0;
}
