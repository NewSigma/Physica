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
#include "Physica/Core/Physics/MD/ForceModel/Ewald/RSpaceEwald.h"

using namespace Physica;

// FIXME: Re-enable this test
constexpr bool Disable = sizeof(int) == 0;
namespace Physica {
    class Test {
        using ValueType = float64;
        using ScalarType = Diff<ValueType, DiffMode::Reverse, 1>;
        using MDCellType = MDCell<ScalarType>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using MassVector = MDCellType::MassVector;
    public:
        static void run() {
            if constexpr (Disable) {
                ScalarType volume = 125;
                const size_t cellSize = 3;
                const auto cell = makeSystem(volume, cellSize);
                const auto& pos = cell.getPos();
                VectorND<ScalarType> charges(cell.getNumParticle(), 1.0);
                auto tail = charges.tail(cell.getNumParticle() / 2);
                tail = float64(-1);
                RSpaceEwald<ScalarType, false> ewald(cell.getLattice(), std::move(charges));
                ewald.potentialV(pos).reverse();
                /* Test press */ {
                    const ValueType press_diff = -volume.grad() / ValueType(cellSize * cellSize * cellSize);
                    const ValueType press = (ewald.virial(pos).trace() / ScalarType(3)).value();
                    if (!scalarNear(press_diff, press, 1E-13))
                        exit(EXIT_FAILURE);
                }
            }
        }
    private:
        static MDCellType makeSystem(ScalarType& volume, size_t cellSize) {
            if constexpr (Disable) {
                constexpr size_t numMolecularUnitCell = 2;
                LatticeMatrix lattice = MDCellType::LatticeMatrix::identity(3);
                const auto latticeConst = cbrt(volume);
                lattice *= latticeConst;

                PositionMatrix pos(numMolecularUnitCell, 3, ValueType(0));
                auto row = pos.row(1);
                row = float64(0.5);
                pos *= latticeConst;

                MassVector massVec(numMolecularUnitCell, 1);
                MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));
                cell.toSuperCell<ExtendCellOption::AtomMajor>({cellSize, cellSize, cellSize});
                return cell;
            }
            return {};
        }
    };
}

int main() {
    Physica::Test::run();
    return 0;
}
