/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Physics/ElectronStructure/CrystalCell.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "WaveFunction.h"
#include "KPointGrid.h"

namespace Physica::Core {
    template<class ScalarType>
    class KSSolver {
        using KPoint = typename KPointGrid::KPoint;
        using Hamilton = DenseSymmMatrix<ScalarType>;
        using KSOrbit = WaveFunction<ScalarType>;

        CrystalCell cell;
        ScalarType cutEnergy;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergy_);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver() = default;
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte);
    private:
        static void fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
        static void fillExternalPot(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
        static void fillHartree(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
        static void fillXC(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
    };

    template<class ScalarType>
    KSSolver<ScalarType>::KSSolver(CrystalCell cell_, ScalarType cutEnergy_)
            : cell(std::move(cell_))
            , cutEnergy(std::move(cutEnergy_)) {}

    template<class ScalarType>
    bool KSSolver::solve(const ScalarType& criteria, size_t maxIte) {
        KSOrbit orbit(cutEnergy, cell.getReciprocal());
        Hamilton hamilton = Hamilton(orbit.getPlainWaveCount());
        KPoint toSolve{0, 0, 0};
        while (true) {

        };
        return true;
    }

    template<class ScalarType>
    void KSSolver::fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit) {
        const size_t order = hamilton.getRow();
        for (size_t i = 0; i < order; ++i)
            hamilton(i, i) += (k + orbit.getBaseFunc(i)).squaredNorm() * ScalarType(0.5);
    }
}
