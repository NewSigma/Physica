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

#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"
#include "Physica/Core/Physics/ElectronStructure/CrystalCell.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "WaveFunction.h"
#include "KPointGrid.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType>
    class KSSolver {
        using KPoint = typename KPointGrid::KPoint;
        using Hamilton = DenseSymmMatrix<ScalarType>;
        using KSOrbit = WaveFunction<ScalarType>;
        using KSOrbits = Utils::Array<KSOrbit>;

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
        [[nodiscard]] static size_t numOrbitToSolve(const CrystalCell& cell);
        static void fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
        static void fillPotential(KPoint k, Hamilton& hamilton, const Grid3D& chargeDensity);
        void sortEigenvalues(Utils::Array<size_t>& indexToSort) const;
        static void updateOrbits(const EigenSolver<Hamilton>& eigenSolver, KSOrbits& orbits);
        static void updateDensity(Grid3D& chargeDensity, const KSOrbits& orbits);
    };

    template<class ScalarType>
    KSSolver<ScalarType>::KSSolver(CrystalCell cell_, ScalarType cutEnergy_)
            : cell(std::move(cell_))
            , cutEnergy(std::move(cutEnergy_)) {}

    template<class ScalarType>
    bool KSSolver<ScalarType>::solve(const ScalarType& criteria, size_t maxIte) {
        const size_t orbitCount = numOrbitToSolve(cell);
        auto orbits = KSOrbits(numOrbitToSolve(cell), KSOrbit(cutEnergy, cell.getReciprocal()));
        auto chargeDensity = Grid3D<ScalarType>();
        KPoint toSolve{0, 0, 0};

        const size_t plainWaveCount = orbit.getPlainWaveCount();
        auto hamilton = Hamilton(plainWaveCount);
        auto eigenSolver = EigenSolver<Hamilton>(plainWaveCount);
        auto sortedEigenvalues = Utils::Array<size_t>(orbitCount);

        size_t iteration = 0;
        while (true) {
            fillKinetic(toSolve, hamilton, orbits[0]); //Any orbit is ok, we need base function only
            fillPotential(toSolve, hamilton, chargeDensity);
            eigenSolver.compute(hamilton, true);
            
            sortEigenvalues(sortedEigenvalues);
            if (++iteration >= matIte)
                break;
            hamilton = ScalarType::Zero();
            updateOrbits(eigenSolver, orbits);
            updateDensity(chargeDensity, orbits);
        };
        return true;
    }

    template<class ScalarType>
    size_t KSSolver<ScalarType>::numOrbitToSolve(const CrystalCell& cell) {
        size_t result = 0;
        for (size_t i = 0; i < cell.getAtomCount(); ++i)
            result += cell.getCharge(i)
        return result;
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit) {
        const size_t order = hamilton.getRow();
        for (size_t i = 0; i < order; ++i)
            hamilton(i, i) += (k + orbit.getBaseFunc(i)).squaredNorm() * ScalarType(0.5);
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::fillPotential(KPoint k, Hamilton& hamilton, const Grid3D& chargeDensity) {
    }
    /**
     * Reference to RHFSolver module
     * TODO: eliminate code duplicate
     */
    template<class ScalarType>
    void KSSolver<ScalarType>::sortEigenvalues(Utils::Array<size_t>& indexToSort) const {
        const auto& eigenvalues = eigenSolver.getEigenvalues();
        auto arrayToSort = toRealVector(eigenvalues);
        const size_t length = indexToSort.getLength();
        for (size_t i = 0; i < length; ++i)
            indexToSort[i] = i;

        for (size_t i = 0; i < length; ++i) {
            size_t indexOfToInsert = i;
            for (size_t j = i + 1; j < length; ++j) {
                if (arrayToSort[indexOfToInsert] > arrayToSort[j])
                    indexOfToInsert = j;
            }
            std::swap(arrayToSort[i], arrayToSort[indexOfToInsert]);
            std::swap(indexToSort[i], indexToSort[indexOfToInsert]);
            assert(eigenvalues[indexToSort[i]].getReal() >= eigenvalues[indexToSort[i == 0 ? 0 : i - 1]].getReal());
        }
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::updateOrbits(const EigenSolver<Hamilton>& eigenSolver, KSOrbits& orbits) {
        const auto eigenVectors = eigenSolver.getEigenvectors();
        const size_t orbitCount = orbits.getLength();
        for (size_t i = 0; i < orbitCount; ++i)
            orbits[i] = toRealVector(eigenVectors.col(i));
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::updateDensity(Grid3D& chargeDensity, const KSOrbits& orbits) {
        auto[dimX, dimY, dimZ] = chargeDensity.getDim();
        for (size_t i = 0; i < dimX; ++i) {
            for (size_t j = 0; j < dimY; ++j) {
                for (size_t k = 0; k < dimZ; ++j) {
                    const auto pos = chargeDensity.dimToPos({i, j, k});
                    auto density = ScalarType::Zero();
                    for (const auto& orbit : orbits)
                        density += orbit(pos).squaredNorm();
                    chargeDensity(i, j, k) = density;
                }
            }
        }
    }
}
