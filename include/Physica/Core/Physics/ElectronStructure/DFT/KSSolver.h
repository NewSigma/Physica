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

#include "Physica/Core/Exception/BadConvergenceException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
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
        using Hamilton = DenseSymmMatrix<ScalarType>;
        using KSOrbit = WaveFunction<ScalarType>;
        using KSOrbits = Utils::Array<KSOrbit>;
        using UnsignedGrid = Grid3D<ScalarType, false>;
        using SignedGrid = Grid3D<ComplexType, true>;

        CrystalCell cell;
        ReciprocalCell repCell;
        ScalarType cutEnergy;
        UnsignedGrid densityGrid;
        UnsignedGrid totalPotGrid;
        SignedGrid externalPotGrid;
        FFT<ScalarType, 3> fftSolver;
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
        [[nodiscard]] size_t electronCount() const;
        [[nodiscard]] size_t numOrbitToSolve() const;
        void initDensity();
        void initExternalPot();
        static void fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit);
        void fillPotential(Hamilton& hamilton, const KSOrbit& orbit);
        static void updateOrbits(const EigenSolver<Hamilton>& eigenSolver, KSOrbits& orbits);
        void updateDensity(const KSOrbits& orbits);
        void updateHartree();
        void updateXCPot();
        [[nodiscard]] Utils::Array<SignedGrid> getStructureFactor();
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType>
    KSSolver<ScalarType>::KSSolver(CrystalCell cell_, ScalarType cutEnergy_, size_t gridDimX, size_t gridDimY, size_t gridDimZ)
            : cell(std::move(cell_))
            , repCell(cell_.reciprocal())
            , cutEnergy(std::move(cutEnergy_))
            , densityGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ)
            , totalPotGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ)
            , fftSolver({gridDimX, gridDimY, gridDimZ}, {ScalarType(cell_.getLattice().row(0).norm()) / ScalarType(gridDimX - 1),
                                                         ScalarType(cell_.getLattice().row(1).norm()) / ScalarType(gridDimY - 1),
                                                         ScalarType(cell_.getLattice().row(2).norm()) / ScalarType(gridDimZ - 1)}) {
        externalPotGrid = SignedGrid::gridFromCutEnergy(cutEnergy, repCell.getLattice());
        initExternalPot();
    }

    template<class ScalarType>
    bool KSSolver<ScalarType>::solve(const ScalarType& criteria, size_t maxIte) {
        auto orbits = KSOrbits(numOrbitToSolve(), KSOrbit(cutEnergy, repCell.getLattice()));
        KPoint toSolve{0, 0, 0};

        const size_t plainWaveCount = orbits[0].getPlainWaveCount();
        const ScalarType inv_volume = reciprocal(cell.getVolume());
        auto hamilton = Hamilton(plainWaveCount);
        auto eigenSolver = EigenSolver<Hamilton>(plainWaveCount);
        UnsignedGrid lastDensity = densityGrid;

        size_t iteration = 0;
        while (true) {
            fillKinetic(toSolve, hamilton, orbits[0]); //Any orbit is ok, we need base function only
            fillPotential(hamilton, orbits[0]);
            eigenSolver.compute(hamilton, true);
            eigenSolver.sort();

            hamilton = ScalarType::Zero();
            updateOrbits(eigenSolver, orbits);
            updateDensity(orbits);

            const bool isDensityConverged = abs(divide((densityGrid.asVector() - lastDensity.asVector()), densityGrid.asVector())).max() < criteria;
            if (isDensityConverged)
                break;

            if (++iteration == maxIte)
                throw BadConvergenceException();
        };
        return true;
    }

    template<class ScalarType>
    size_t KSSolver<ScalarType>::electronCount() const {
        size_t result = 0;
        for (size_t i = 0; i < cell.getAtomCount(); ++i)
            result += getCharge(cell.getAtomicNumber(i));
        return result;
    }

    template<class ScalarType>
    size_t KSSolver<ScalarType>::numOrbitToSolve() const {
        return (electronCount() + 1) / 2;
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::initDensity() {
        densityGrid.asVector() = ScalarType(electronCount(cell)) / cell.getVolume();
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::initExternalPot() {
        const Utils::Array<SignedGrid> all_factors = getStructureFactor();
        const ScalarType factor1 = ScalarType(-4 * M_PI) / cell.getVolume();
        const std::unordered_set<uint16_t> species = cell.getSpecies();

        externalPotGrid.asVector() = ScalarType::Zero();
        const size_t gridSize = externalPotGrid.getSize();

        size_t j = 0;
        for (uint16_t element : species) {
            const SignedGrid& factors = all_factors[j];
            for (size_t i = 0; i < gridSize; ++i)
                externalPotGrid[i] += factor1 * getCharge(element) * factors[i] / factors.indexToPos(i).squaredNorm();
            ++j;
        }
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::fillKinetic(KPoint k, Hamilton& hamilton, const KSOrbit& orbit) {
        const size_t order = hamilton.getRow();
        for (size_t i = 0; i < order; ++i)
            hamilton(i, i) += ScalarType((k + orbit.getWaveVector(orbit.indexToDim(i))).squaredNorm()) * ScalarType(0.5);
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::fillPotential(Hamilton& hamilton, const KSOrbit& orbit) {
        totalPotGrid.asVector() = ScalarType::Zero();
        updateHartree();
        updateXCPot();
        fftSolver.transform(totalPotGrid.asVector());

        const size_t order = hamilton.getRow();
        for (size_t i = 0; i < order; ++i) {
            const Vector<ScalarType, 3> k1 = orbit.getWaveVector(orbit.indexToDim(i));
            for (size_t j = i; j < order; ++j) {
                const Vector<ScalarType, 3> k2 = orbit.getWaveVector(orbit.indexToDim(j));
                hamilton(i, j) += fftSolver.getFreqIntense(Vector<ScalarType, 3>((k1 - k2) / ScalarType(2 * M_PI))).norm(); //Not norm
            }
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
    void KSSolver<ScalarType>::updateDensity(const KSOrbits& orbits) {
        auto[dimX, dimY, dimZ] = densityGrid.getDim();
        for (size_t i = 0; i < dimX; ++i) {
            for (size_t j = 0; j < dimY; ++j) {
                for (size_t k = 0; k < dimZ; ++j) {
                    const auto pos = densityGrid.dimToPos({i, j, k});
                    auto density = ScalarType::Zero();
                    for (size_t index = 0; index < orbits.getLength(); ++index)
                        density += orbits[index](pos).squaredNorm(); //We need multiply the occupation number
                    densityGrid(i, j, k) = density;
                }
            }
        }
    }

    template<class ScalarType>
    void KSSolver<ScalarType>::updateHartree() {
        const size_t gridSize = totalPotGrid.getSize();
        for (size_t i = 0; i < gridSize; ++i)
            totalPotGrid[i] += Ewald<ScalarType>::potHartree(totalPotGrid.indexToPos(i), densityGrid, repCell);
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
        totalPotGrid.asVector() += pow(densityGrid.asVector(), ScalarType(1.0 / 3)) * ScalarType(exchange_factor);
        totalPotGrid.asVector() += ScalarType(correlation_factor1) * ln(ScalarType::One() + ScalarType(correlation_factor2) * pow(densityGrid.asVector(), ScalarType(1.0 / 3))); //Reference [1]
    }

    template<class ScalarType>
    typename Utils::Array<typename KSSolver<ScalarType>::SignedGrid> KSSolver<ScalarType>::getStructureFactor() {
        const std::unordered_set<uint16_t> species = cell.getSpecies();
        const auto& lattice = repCell.getLattice();
        auto all_factors = Utils::Array<SignedGrid>(species.size(), SignedGrid::gridFromCutEnergy(cutEnergy, lattice));
        const size_t factors_size = all_factors[0].getSize();
        const size_t atomCount = cell.getAtomCount();

        Vector<ScalarType, 3> g;
        size_t j = 0;
        for (int16_t element : species) {
            SignedGrid& factors = all_factors[j];
            for (size_t i = 0; i < factors_size; ++i) {
                g = factors.indexToPos(i);
                auto temp = ComplexType::Zero();
                for (size_t ion = 0; ion < atomCount; ++i) {
                    if (cell.getAtomicNumber(ion) == element) {
                        auto r = cell.getPos().row(i);
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
