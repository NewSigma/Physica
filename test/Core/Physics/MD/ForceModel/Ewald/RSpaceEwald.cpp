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
#include "Physica/Core/Physics/MD/ForceModel/Ewald/RSpaceEwald.h"
#include "Test.h"

using namespace Physica;

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
    public:
        static void run() {
            if constexpr (Disable) {
                dfloat volume = 125;
                const size_t cellSize = 3;
                const auto cell = makeSystem(volume, cellSize);
                const auto& pos = cell.getPos();
                VectorND<dfloat> charges(cell.getNumParticle(), 1.0);
                auto tail = charges.tail(cell.getNumParticle() / 2);
                tail = T(-1);
                RSpaceEwald<dfloat, false> ewald(cell.getLattice(), std::move(charges));
                ewald.potentialV(pos).reverse();
                /* Test press */ {
                    const T press_diff = -volume.grad() / T(cellSize * cellSize * cellSize);
                    const T press = (ewald.virial(pos).trace() / dfloat(3)).value();
                    expect(scalarNear(press_diff, press, 1E-13));
                }
            }
        }
    private:
        static MDCellType makeSystem(dfloat& volume, size_t cellSize) {
            if constexpr (Disable) {
                constexpr size_t numMolecularUnitCell = 2;
                LatticeMatrix lattice = MDCellType::LatticeMatrix::identity(3);
                const auto latticeConst = cbrt(volume);
                lattice *= latticeConst.value();

                PositionMatrix pos(numMolecularUnitCell, 3, T(0));
                auto row = pos.row(1);
                row = T(0.5);
                pos *= latticeConst.value();

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
