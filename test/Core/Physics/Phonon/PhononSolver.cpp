/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/Phonon/PhononSolver.h"
#include "Physica/Core/Physics/Phonon/PhononSolverImpl/FCProjector.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using SolverType = PhononSolver<T>;
using GridType = SolverType::RSpaceFCGrid;

namespace {
    void projector() {
        const Index3D superSize{2, 2, 2};
        const size_t numDOF = 3;
        const FCProjector<T> projector(superSize, numDOF);

        auto grid = GridType::zeros(superSize[0], superSize[1], superSize[2], numDOF, numDOF);
        grid[0UZ, 0UZ, 0UZ, 0UZ, 0UZ] = T(1);
        grid[1UZ, 0UZ, 1UZ, 2UZ, 1UZ] = T(-2);

        const auto vec = projector.toVector(grid);
        expect(vec.getLength() == projector.getNumForceConsts());

        auto result = GridType::zeros(superSize[0], superSize[1], superSize[2], numDOF, numDOF);
        projector.toGrid(vec, result);
        expect(result == grid);
    }

    void toKSpace() {
        using MDCellType = SolverType::MDCellType;
        const MDCellType unitCell(1, MDCellType::MassVector(1, T(1)));
        const SolverType solver(unitCell, Index3D{2, 2, 2});

        const size_t dof = solver.getUnitCellDOF();
        auto grid = GridType::zeros(2UZ, 2UZ, 2UZ, dof, dof);
        grid[0UZ, 0UZ, 0UZ, 0UZ, 0UZ] = T(1);
        grid[0UZ, 0UZ, 0UZ, 1UZ, 2UZ] = T(-2);

        const auto kSpace = solver.toKSpace(grid);
        expect(kSpace.dim(0) == 2UZ && kSpace.dim(1) == 2UZ && kSpace.dim(2) == 2UZ);
        expect(kSpace.dim(3) == dof && kSpace.dim(4) == dof);

        // A delta in real space transforms to a constant in k space
        for (size_t i = 0; i < 2UZ; ++i)
            for (size_t j = 0; j < 2UZ; ++j)
                for (size_t k = 0; k < 2UZ; ++k) {
                    expect(scalarNear(kSpace[i, j, k, 0UZ, 0UZ].real(), T(1), 1E-12));
                    expect(abs(kSpace[i, j, k, 0UZ, 0UZ].imag()) < 1E-12);
                    expect(scalarNear(kSpace[i, j, k, 1UZ, 2UZ].real(), T(-2), 1E-12));
                    expect(abs(kSpace[i, j, k, 1UZ, 2UZ].imag()) < 1E-12);
                    expect(abs(kSpace[i, j, k, 0UZ, 1UZ]) < 1E-12);
                }
    }
}

int main() {
    projector();
    toKSpace();
    return 0;
}
