/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Physica/Core/Exception/BadConvergenceException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/ElectronStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronStructure/ReciprocalCell.h"
#include "Ewald.h"
#include "WaveFunction.h"
#include "KPointGrid.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType>
    class KSSolver {
        using ComplexType = ComplexScalar<ScalarType>;
        using KPoint = typename KPointGrid::KPoint;
        using Hamilton = DenseHermiteMatrix<ComplexType>;
        using MatrixType = DenseMatrix<ComplexType>;
        using KSOrbit = WaveFunction<ScalarType>;
        using KSOrbits = Utils::Array<KSOrbit>;
        using UncenteredGrid = Grid3D<ScalarType, false>;
        using CenteredGrid = Grid3D<ComplexType, true>;

        CrystalCell cell;
        ReciprocalCell repCell;
        ScalarType cutEnergy;
        UncenteredGrid densityGrid;
        UncenteredGrid xcPotGrid;
        CenteredGrid externalPotGrid;
        FFT<ScalarType, 3> fftSolver;
        size_t iteration;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergy_, size_t gridDimX_, size_t gridDimY_, size_t gridDimZ_);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver() = default;
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte);
    private:
        [[nodiscard]] size_t numOrbitToSolve() const;
        [[nodiscard]] ScalarType occupacy(size_t orbitIndex) const;
        void initDensity();
        void initExternalPot();
        static void fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
        void fillPotential(Hamilton& hamilton, const KSOrbit& orbit);
        static void updateOrbits(const EigenSolver<MatrixType>& eigenSolver, KSOrbits& orbits);
        void updateDensity(const KSOrbits& orbits);
        void updateXCPot();
        [[nodiscard]] Utils::Array<CenteredGrid> getStructureFactor(ScalarType factorCutoff);
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType>
    KSSolver<ScalarType>::KSSolver(CrystalCell cell_, ScalarType cutEnergy_, size_t gridDimX, size_t gridDimY, size_t gridDimZ)
            : cell(std::move(cell_))
            , repCell(cell_.reciprocal())
            , cutEnergy(std::move(cutEnergy_))
            , densityGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ)
            , xcPotGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ)
            , fftSolver({gridDimX, gridDimY, gridDimZ}, {ScalarType(cell_.getLattice().row(0).norm()) / ScalarType(gridDimX - 1),
                                                         ScalarType(cell_.getLattice().row(1).norm()) / ScalarType(gridDimY - 1),
                                                         ScalarType(cell_.getLattice().row(2).norm()) / ScalarType(gridDimZ - 1)})
            , iteration(0) {
        initExternalPot();
    }

    template<class ScalarType>
    bool KSSolver<ScalarType>::solve(const ScalarType& criteria, size_t maxIte) {
        auto orbits = KSOrbits(numOrbitToSolve(), KSOrbit(cutEnergy, repCell.getLattice()));
        KPoint toSolve{0, 0, 0};

        const size_t plainWaveCount = orbits[0].getPlainWaveCount();
        const ScalarType inv_volume = reciprocal(cell.getVolume());
        auto hamilton = Hamilton(plainWaveCount);
        auto eigenSolver = EigenSolver<MatrixType>(plainWaveCount);
        UncenteredGrid lastDensity = densityGrid;

        iteration = 0;
        while (true) {
            fillKinetic(toSolve, hamilton, orbits[0]); //Any orbit is ok, we need base function only
            fillPotential(hamilton, orbits[0]);
            eigenSolver.compute(hamilton, true);
            eigenSolver.sort();

            hamilton = ScalarType::Zero();
            updateOrbits(eigenSolver, orbits);
            updateDensity(orbits);

            const ScalarType densityChange = abs(divide((densityGrid.asVector() - lastDensity.asVector()), densityGrid.asVector())).max();
            lastDensity = densityGrid;
            const bool isDensityConverged = densityChange < criteria;
            if (isDensityConverged)
                break;

            if (++iteration == maxIte)
                throw BadConvergenceException();
        };
        return true;
    }

    template<class ScalarType>
    size_t KSSolver<ScalarType>::numOrbitToSolve() const {
        return (cell.getElectronCount() + 1) / 2;
    }

    template<class ScalarType>
    ScalarType KSSolver<ScalarType>::occupacy(size_t orbitIndex) const {
        if (orbitIndex == numOrbitToSolve() - 1)
            return ScalarType(2 - (cell.getElectronCount() % 2U == 0U));
        else
            return ScalarType::Two();
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::initDensity() {
        densityGrid.asVector() = ScalarType(cell.getElectronCount()) / cell.getVolume();
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::initExternalPot() {
        const ScalarType factorCutoff = cutEnergy * 8;
        externalPotGrid = CenteredGrid::gridFromCutEnergy(factorCutoff, repCell.getLattice());
        const Utils::Array<CenteredGrid> all_factors = getStructureFactor(factorCutoff);
        const ScalarType factor1 = ScalarType(-4 * M_PI) / cell.getVolume();
        const std::unordered_set<uint16_t> species = cell.getSpecies();

        externalPotGrid.asVector() = ScalarType::Zero();
        const size_t gridSize = externalPotGrid.getSize();

        size_t j = 0;
        for (uint16_t element : species) {
            const CenteredGrid& factors = all_factors[j];
            for (size_t i = 0; i < gridSize; ++i)
                externalPotGrid[i] += factor1 * getCharge(element) * factors[i] / factors.indexToPos(i).squaredNorm();
            ++j;
        }
        externalPotGrid(0, 0, 0) = ComplexType::Zero();
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit) {
        const size_t order = hamilton.getRow();
        for (size_t i = 0; i < order; ++i)
            hamilton(i, i) += ScalarType((k + orbit.getWaveVector(orbit.indexToDim(i))).squaredNorm()) * ScalarType(0.5);
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::fillPotential(Hamilton& hamilton, const KSOrbit& orbit) {
        using Dim = typename CenteredGrid::Dim;
        xcPotGrid.asVector() = ScalarType::Zero();
        updateXCPot();
        fftSolver.transform(xcPotGrid.asVector());

        const ScalarType factor1 = reciprocal(ScalarType(2 * M_PI));
        const size_t order = hamilton.getRow();
        for (size_t i = 0; i < order; ++i) {
            const Dim dim1 = orbit.indexToDim(i);
            auto[x1, y1, z1] = dim1;
            const Vector<ScalarType, 3> k1 = orbit.getWaveVector(dim1);
            for (size_t j = i; j < order; ++j) {
                const Dim dim2 = orbit.indexToDim(j);
                auto[x2, y2, z2] = dim2;
                const Vector<ScalarType, 3> k2 = orbit.getWaveVector(dim2);
                hamilton(i, j) += fftSolver.getFreqIntense(Vector<ScalarType, 3>((k1 - k2) * factor1)) + externalPotGrid(x1 - x2, y1 - y2, z1 - z2);
            }
        }

        fftSolver.transform(densityGrid.asVector());
        const ScalarType factor = ScalarType(4 * M_PI) / cell.getVolume();
        for (size_t i = 0; i < order; ++i) {
            const Vector<ScalarType, 3> k1 = orbit.getWaveVector(orbit.indexToDim(i));
            for (size_t j = i + 1; j < order; ++j) {
                const Vector<ScalarType, 3> k2 = orbit.getWaveVector(orbit.indexToDim(j));
                const Vector<ScalarType, 3> k = k1 - k2;
                hamilton(i, j) += fftSolver.getFreqIntense(Vector<ScalarType, 3>(k * factor1)) * factor / k.squaredNorm();
            }
        }
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::updateOrbits(const EigenSolver<MatrixType>& eigenSolver, KSOrbits& orbits) {
        const auto& eigenVectors = eigenSolver.getRawEigenvectors();
        const size_t orbitCount = orbits.getLength();
        for (size_t i = 0; i < orbitCount; ++i)
            orbits[i] = eigenVectors.col(i);
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::updateDensity(const KSOrbits& orbits) {
        auto[dimX, dimY, dimZ] = densityGrid.getDim();
        for (size_t i = 0; i < dimX; ++i) {
            for (size_t j = 0; j < dimY; ++j) {
                for (size_t k = 0; k < dimZ; ++k) {
                    const auto pos = densityGrid.dimToPos({i, j, k});
                    auto density = ScalarType::Zero();
                    for (size_t index = 0; index < orbits.getLength(); ++index)
                        density += orbits[index](pos).squaredNorm() * occupacy(index);
                    densityGrid(i, j, k) = density;
                }
            }
        }
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:479
     */
    template<class ScalarType>
    void KSSolver<ScalarType>::updateXCPot() {
        constexpr double exchange_factor = -0.98474502184269654;
        constexpr double correlation_factor1 = -0.045 / 2;
        constexpr double correlation_factor2 = 33.851831034345862;
        xcPotGrid.asVector() += pow(densityGrid.asVector(), ScalarType(1.0 / 3)) * ScalarType(exchange_factor);
        xcPotGrid.asVector() += ScalarType(correlation_factor1) * ln(ScalarType::One() + ScalarType(correlation_factor2) * pow(densityGrid.asVector(), ScalarType(1.0 / 3))); //Reference [1]
    }

    template<class ScalarType>
    typename Utils::Array<typename KSSolver<ScalarType>::CenteredGrid> KSSolver<ScalarType>::getStructureFactor(ScalarType factorCutoff) {
        const std::unordered_set<uint16_t> species = cell.getSpecies();
        const auto& lattice = repCell.getLattice();
        auto all_factors = Utils::Array<CenteredGrid>(species.size(), CenteredGrid::gridFromCutEnergy(factorCutoff, lattice));
        const size_t factors_size = all_factors[0].getSize();
        const size_t atomCount = cell.getAtomCount();

        Vector<ScalarType, 3> g;
        size_t j = 0;
        for (uint16_t element : species) {
            CenteredGrid& factors = all_factors[j];
            for (size_t i = 0; i < factors_size; ++i) {
                g = factors.indexToPos(i);
                auto temp = ComplexType::Zero();
                for (size_t ion = 0; ion < atomCount; ++ion) {
                    if (cell.getAtomicNumber(ion) == element) { //We can use searching table method
                        auto r = cell.getPos().row(ion);
                        const ScalarType phase = g * r;
                        temp += ComplexType(cos(phase), sin(phase));
                    }
                }
                factors.asVector()[i] = temp;
            }
            ++j;
        }
        return all_factors;
    }
}
