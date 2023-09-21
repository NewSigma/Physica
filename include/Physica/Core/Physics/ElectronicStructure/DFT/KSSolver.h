/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald.h"
#include "Physica/Utils/TestHelper.h"
#include "Basis/PlainWaveBasis.h"
#include "BandGrid.h"
#include "ChargeMixer.h"
#include "KSHamilton.h"

namespace Physica::Core {
    template<class ScalarType, class XCProvider>
    class KSSolver {
    public:
        constexpr static bool IsSpinPolarized = XCProvider::IsSpinPolarized;
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using RepCellType = ReciprocalCell<typename CrystalCell::ScalarType>;
        using BandType = BandGrid<ScalarType, IsSpinPolarized>;
        using HamiltonType = KSHamilton<ScalarType, IsSpinPolarized>;
        using EigenSolverType = SpinPair<JacobiDavidson<ComplexType>, IsSpinPolarized>;
        using DensityType = DensityGrid<ScalarType, IsSpinPolarized>;
        using BasisType = typename DensityType::BasisType;
        using KSOrbitArray = typename DensityType::KSOrbitArray;
        using PotType = SpinPair<RSpaceGrid<ScalarType>, IsSpinPolarized>;
        using Index3D = typename GridBase::Index3D;

        using LatticeMatrix = typename HamiltonType::LatticeMatrix;
        using FFT3D = FFT<ComplexType, 3>;
    protected:
        size_t iteration;

        KSOrbitArray orbits;
        HamiltonType hamiltonH;
        BandType band;

        EigenSolverType eigSolver;
        DensityType density;
        ChargeMixer<ScalarType, IsSpinPolarized> chargeMixer;
        PotType xcPot;
        XCProvider xcProvider;

        FFT3D fft;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergyPsi_, ScalarType cutEnergyRho_, BandType band_, size_t numBand);
        KSSolver(const KSSolver&) = default;
        KSSolver(KSSolver&&) noexcept = default;
        ~KSSolver() = default;
        /* Operators */
        KSSolver& operator=(KSSolver obj) noexcept;
        /* Operations */
        bool solve(ScalarType criteria, size_t maxIte);

        template<class RandomGenerator> void initWaveFunc(const SpinPair<BasisType, IsSpinPolarized>& initial, RandomGenerator& gen);
        void swap(KSSolver& obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return hamiltonH.getRepLattice(); }
        [[nodiscard]] const LatticeMatrix& getRepLattice() const noexcept { return hamiltonH.getRepLattice(); }
        [[nodiscard]] size_t getNumElectron() const noexcept { return hamiltonH.getCrystalCell().getElectronCount(); }
        [[nodiscard]] size_t getNumPlainWave() const noexcept { return orbits[0][SpinState::Up].getNumPlainWave(); }
        [[nodiscard]] size_t getNumBand() const noexcept { return orbits.getLength(); }
        [[nodiscard]] size_t getNumFilledBand() const noexcept { return (getNumElectron() + 1) / 2; }
        [[nodiscard]] KSOrbitArray& getOrbits() noexcept { return orbits; }
        [[nodiscard]] const KSOrbitArray& getOrbits() const noexcept { return orbits; }
        [[nodiscard]] const BandType& getBand() const noexcept { return band; }
        [[nodiscard]] DensityType& getDensity() noexcept { return density; }
        [[nodiscard]] const DensityType& getDensity() const noexcept { return density; }
    protected:
        /* Operations */
        void calcDensity(const BasisType& orbit);
        void updateOrbits();
        void updateDensity();

        friend class Physica::Test;
    };

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::KSSolver(
                CrystalCell cell_,
                ScalarType cutEnergyPsi_,
                ScalarType cutEnergyRho_,
                BandType band_,
                size_t numBand)
            : iteration(0)
            , hamiltonH(std::move(cell_), cutEnergyPsi_, cutEnergyRho_)
            , band(std::move(band_)) {
        assert(numBand >= getNumFilledBand() && "[Error]: Number of band is not enough to hold all electrons");
        orbits.resize(numBand, hamiltonH.getCutEnergyPsi(), getRepLattice());
        /* Eigensolver */ {
            const size_t planeWaveCount = hamiltonH.getNumBasis();
            eigSolver = EigenSolverType();
            {
                eigSolver[SpinState::Up].resize(planeWaveCount, numBand);
            }
            if constexpr (IsSpinPolarized) {
                eigSolver[SpinState::Down].resize(planeWaveCount, numBand);
            }
        }
        const auto rSpaceSize = hamiltonH.getFFTRSpaceSize();
        density.resize(rSpaceSize);

        fft = FFT3D(rSpaceSize, {1, 1, 1}, PlanFlag::Estimate);

        chargeMixer = ChargeMixer<ScalarType, IsSpinPolarized>(getRepLattice(), rSpaceSize);
        xcPot = PotType(rSpaceSize);
        xcProvider = XCProvider(getNumPlainWave());
    }

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>& KSSolver<ScalarType, XCProvider>::operator=(KSSolver obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class XCProvider>
    bool KSSolver<ScalarType, XCProvider>::solve(ScalarType criteria, size_t maxIte) {
        for (auto& kPoint : band.getKPointGrid()) {
            iteration = 0;
            while (true) {
                const Vector3D k = kPoint.getPos();
                hamiltonH.makeHamiltonWithXC(xcPot[SpinState::Up], density, kPoint.getPos());
                /* Solve band */ {
                    {
                        const auto orbit = orbits[0][SpinState::Up].getCoeffGrid().flatten();
                        eigSolver[SpinState::Up].compute(hamiltonH[SpinState::Up], orbit);
                        eigSolver[SpinState::Up].sort();
                    }
                    if constexpr (IsSpinPolarized) {
                        const auto orbit = orbits[0][SpinState::Down].getCoeffGrid().flatten();
                        eigSolver[SpinState::Down].compute(hamiltonH[SpinState::Down], orbit);
                        eigSolver[SpinState::Down].sort();
                    }
                }
                updateOrbits();
                updateDensity();

                if (iteration != 0) {
                    const auto& delta_rho = chargeMixer.getResidule().getTotalDensity().flatten();
                    const ScalarType error = delta_rho.squaredNorm();
                    const bool isConverged = error < criteria;
                    std::cout << error << std::endl;
                    if (isConverged)
                        break;
                }
                if (iteration == maxIte)
                    throw BadConvergenceException("Exceed max iteration of KSSolver");
                chargeMixer.mix(iteration, density);
                iteration += 1;
            };
            kPoint.setBandEnergy(SpinState::Up, toRealVector(eigSolver[SpinState::Up].getEigenvalues()));
            if constexpr (IsSpinPolarized)
                kPoint.setBandEnergy(SpinState::Down, toRealVector(eigSolver[SpinState::Down].getEigenvalues()));
        }
        return true;
    }

    template<class ScalarType, class XCProvider>
    template<class RandomGenerator>
    void KSSolver<ScalarType, XCProvider>::initWaveFunc(const SpinPair<BasisType, IsSpinPolarized>& initial, RandomGenerator& gen) {
        const ssize_t dimX = std::min(initial[SpinState::Up].getDimX(), orbits[0][SpinState::Up].getDimX());
        const ssize_t dimY = std::min(initial[SpinState::Up].getDimY(), orbits[0][SpinState::Up].getDimY());
        const ssize_t dimZ = std::min(initial[SpinState::Up].getDimZ(), orbits[0][SpinState::Up].getDimZ());
        {
            auto& waveUp = orbits[0][SpinState::Up];
            waveUp.flatten().random_normal(gen);
            for (ssize_t x = -dimX; x <= dimX; ++x)
                for (ssize_t y = -dimY; y <= dimY; ++y)
                    for (ssize_t z = 0; z <= dimZ; ++z)
                        waveUp(x, y, z) = initial[SpinState::Up](x, y, z);
        }
        if constexpr (IsSpinPolarized) {
            auto& waveDown = orbits[0][SpinState::Down];
            waveDown.flatten().random_normal(gen);
            for (ssize_t x = -dimX; x <= dimX; ++x)
                for (ssize_t y = -dimY; y <= dimY; ++y)
                    for (ssize_t z = 0; z <= dimZ; ++z)
                        waveDown(x, y, z) = initial[SpinState::Down](x, y, z);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::swap(KSSolver& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(iteration, obj.iteration);

        orbits.swap(obj.orbits);
        hamiltonH.swap(obj.hamiltonH);
        band.swap(obj.band);

        eigSolver.swap(obj.eigSolver);
        density.swap(obj.density);
        chargeMixer.swap(obj.chargeMixer);
        xcPot.swap(obj.xcPot);
        xcProvider.swap(obj.xcProvider);

        fft.swap(obj.fft);
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::calcDensity(const BasisType& orbit) {
        const auto& kSpacePsi = orbit.getCoeffGrid();
        const size_t psiDimX = kSpacePsi.getDimX() / 2;
        const size_t psiDimY = kSpacePsi.getDimY() / 2;
        const size_t psiDimZ = kSpacePsi.getDimZ() / 2;

        auto& kSpaceRho = fft.getKSpace();
        kSpaceRho = ScalarType(0);
        const size_t rhoDimX = kSpaceRho.getDimX();
        const size_t rhoDimY = kSpaceRho.getDimY();
        const size_t rhoDimZ = kSpaceRho.getDimZ();
        {
            auto corner = kSpaceRho.leftFrontBottomCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
            corner = kSpacePsi.leftFrontBottomCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.rightFrontBottomCorner({rhoDimX - psiDimX, psiDimY + 1, psiDimZ + 1});
            corner = kSpacePsi.rightFrontBottomCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.leftBackBottomCorner({psiDimX + 1, rhoDimY - psiDimY, psiDimZ + 1});
            corner = kSpacePsi.leftBackBottomCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.rightBackBottomCorner({rhoDimX - psiDimX, rhoDimY - psiDimY, psiDimZ + 1});
            corner = kSpacePsi.rightBackBottomCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.leftFrontTopCorner({psiDimX + 1, psiDimY + 1, rhoDimZ - psiDimZ});
            corner = kSpacePsi.leftFrontTopCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.rightFrontTopCorner({rhoDimX - psiDimX, psiDimY + 1, rhoDimZ - psiDimZ});
            corner = kSpacePsi.rightFrontTopCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.leftBackTopCorner({psiDimX + 1, rhoDimY - psiDimY, rhoDimZ - psiDimZ});
            corner = kSpacePsi.leftBackTopCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        {
            auto corner = kSpaceRho.rightBackTopCorner({rhoDimX - psiDimX, rhoDimY - psiDimY, rhoDimZ - psiDimZ});
            corner = kSpacePsi.rightBackTopCorner({psiDimX + 1, psiDimY + 1, psiDimZ + 1});
        }
        fft.invTransform();
        const size_t numCellOld = kSpacePsi.getSize();
        const size_t numCellNew = kSpaceRho.getSize();
        const ScalarType numElectronScale = sqrt(ScalarType(numCellNew) / ScalarType(numCellOld));
        fft.getRSpace() *= numElectronScale;
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateOrbits() {
        const auto& eigSolverUp = eigSolver[SpinState::Up];
        const auto& eigSolverDown = eigSolver[SpinState::Down];
        for (size_t band = 0; band < getNumFilledBand(); ++band) {
            auto& orbit = orbits[band];
            orbit[SpinState::Up] = eigSolverUp.getEigenvectors().col(band);
            if constexpr (IsSpinPolarized)
                orbit[SpinState::Down] = eigSolverDown.getEigenvectors().col(band);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateDensity() {
        const size_t numFilledBand = getNumFilledBand();
        const int numElectronLastBand = getNumElectron() % 2 == 0 ? 2 : 1;
        chargeMixer.fetchInputDensity(density);
        density.getTotalDensity() = ScalarType(0);
        for (size_t band = 0; band < numFilledBand; ++band) {
            calcDensity(orbits[band][SpinState::Up]);
            const bool isLastBand = band == numFilledBand - 1;
            const int numFillElectron = isLastBand ? numElectronLastBand : 2;
            density.getTotalDensity() += toRealGrid(fft.getRSpace()) * ScalarType(numFillElectron);
        }
        chargeMixer.mix(density);
    }
}
