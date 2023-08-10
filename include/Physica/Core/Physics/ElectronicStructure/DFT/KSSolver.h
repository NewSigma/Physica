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
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "BandGrid.h"
#include "PWBaseWave.h"
#include "Ewald.h"
#include "ChargeMixer.h"

namespace Physica::Core {
    template<class ScalarType, class XCProvider>
    class KSSolver {
    public:
        constexpr static bool isSpinPolarized = XCProvider::isSpinPolarized;
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using RepCellType = ReciprocalCell<typename CrystalCell::ScalarType>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using KSOrbit = PWBaseWave<ScalarType>;
        using BandType = BandGrid<ScalarType, isSpinPolarized>;
        using KSOrbits = SpinPair<Utils::Array<KSOrbit>, isSpinPolarized>;
        using Hamilton = SpinPair<HermiteMatrix, isSpinPolarized>;
        using EigenSolverType = SpinPair<JacobiDavidson<ComplexType>, isSpinPolarized>;
        using DensityType = DensityGrid<ScalarType, isSpinPolarized>;
        using PotType = SpinPair<RSpaceGrid<ScalarType>, isSpinPolarized>;
        using FFT3D = FFT<ScalarType, 3>;
        using FFTxc = SpinPair<FFT3D, isSpinPolarized>;
    protected:
        CrystalCell cell;
        RepCellType repCell;
        ScalarType cutEnergy;
        size_t iteration;

        KSOrbits orbits;
        Hamilton h;
        BandType band;

        EigenSolverType eigSolver;
        DensityType density;
        ChargeMixer<ScalarType, isSpinPolarized> chargeMixer;
        KSpaceGrid<ComplexType> externalPot;
        PotType xcPot;
        XCProvider xcProvider;
        FFT3D fft;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergy_, ScalarType fftDimFactor, BandType band_);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver() = default;
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(ScalarType criteria, size_t maxIte);
        void initDensity();
        void initDensity(const DensityType& rho);
        template<class RandomGenerator> void initWaveFunc(RandomGenerator& gen);
        template<class RandomGenerator> void initWaveFunc(const SpinPair<KSOrbit, isSpinPolarized>& initial, RandomGenerator& gen);
        /* Getters */
        [[nodiscard]] size_t getNumPlaneWave() const noexcept { return orbits[SpinState::Up][0].getNumPlaneWave(); }
        [[nodiscard]] const KSOrbits& getOrbits() const noexcept { return orbits; }
        [[nodiscard]] const BandType& getBand() const noexcept { return band; }
        [[nodiscard]] const DensityType& getDensityGrid() const noexcept { return density; }
        [[nodiscard]] typename RSpaceGrid<ScalarType>::Index3D getFFTDim() const noexcept { return xcPot[SpinState::Up].getDim(); }
    protected:
        /* Operations */
        typename RSpaceGrid<ScalarType>::Index3D makeFFTDim(ScalarType fftDimFactor);
        void initExternalPot();
        void assembleH(Vector3D k);
        void fillPotential();
        void updateOrbits();
        void updateDensity(Vector3D k);
        /* Getters */
        [[nodiscard]] size_t numOrbitToSolve() const { return (cell.getElectronCount() + 1) / 2; }
        [[nodiscard]] size_t getDensitySize() const noexcept { return xcProvider.getBufferSize(); }
        /* Static members */
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::KSSolver(CrystalCell cell_, ScalarType cutEnergy_, ScalarType fftDimFactor, BandType band_)
            : cell(std::move(cell_))
            , cutEnergy(cutEnergy_)
            , iteration(0)
            , band(std::move(band_)) {
        repCell = cell.reciprocal();
        /* Allocate orbits */ {
            const size_t electronCount = cell.getElectronCount();
            orbits[SpinState::Up] = Utils::Array<KSOrbit>((electronCount + 1) / 2, cutEnergy, repCell.getLattice());
            if constexpr (isSpinPolarized)
                orbits[SpinState::Down] = Utils::Array<KSOrbit>(electronCount / 2, cutEnergy, repCell.getLattice());
        }
        /* Eigensolver */ {
            const size_t planeWaveCount = getNumPlaneWave();
            const size_t numBand = band.getNumBand();
            h = Hamilton(planeWaveCount);
            eigSolver = EigenSolverType();
            {
                eigSolver[SpinState::Up].resize(planeWaveCount, numBand);
            }
            if constexpr (isSpinPolarized) {
                eigSolver[SpinState::Down].resize(planeWaveCount, numBand);
            }
        }
        /* Mesh grids */ {
            const auto dim = makeFFTDim(fftDimFactor);
            xcProvider = XCProvider(dim[0] * dim[1] * dim[2]);
            xcPot = PotType(dim);
            density.resize(dim[0], dim[1], dim[2]);
            chargeMixer = ChargeMixer<ScalarType, isSpinPolarized>(repCell, dim);
            /* Allocate fft */ {
                Utils::Array<ScalarType, 3> fftDeltaTs{};
                for (int i = 0; i < 3; ++i)
                    fftDeltaTs[i] = ScalarType(2 * M_PI) / repCell.getLattice().row(i).norm();
                fft = FFT3D(dim, fftDeltaTs, PlanFlag::Estimate);
            }
        }
        initExternalPot();
    }

    template<class ScalarType, class XCProvider>
    bool KSSolver<ScalarType, XCProvider>::solve(ScalarType criteria, size_t maxIte) {
        for (auto& kPoint : band.getKPointGrid()) {
            iteration = 0;
            while (true) {
                const Vector3D k = kPoint.getPos();
                assembleH(k);
                /* Solve band */ {
                    {
                        eigSolver[SpinState::Up].compute(h[SpinState::Up], orbits[SpinState::Up][0].flatten());
                        eigSolver[SpinState::Up].sort();
                    }
                    if constexpr (isSpinPolarized) {
                        eigSolver[SpinState::Down].compute(h[SpinState::Down], orbits[SpinState::Down][0].flatten());
                        eigSolver[SpinState::Down].sort();
                    }
                }
                updateOrbits();
                updateDensity(k);

                if (iteration != 0) {
                    const auto& delta_rho = chargeMixer.getResidule()[SpinState::Up].flatten();
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
            if constexpr (isSpinPolarized)
                kPoint.setBandEnergy(SpinState::Down, toRealVector(eigSolver[SpinState::Down].getEigenvalues()));
        }
        return true;
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initDensity() {
        const ScalarType averageDensity = ScalarType(cell.getElectronCount()) / cell.getVolume();
        {
            auto rho = density[SpinState::Up].flatten();
            rho = averageDensity;
        }
        if constexpr (isSpinPolarized) {
            auto zeta = density[SpinState::Down].flatten();
            zeta = ScalarType(0);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initDensity(const DensityType& rho) {
        density.fit(rho);
    }

    template<class ScalarType, class XCProvider>
    template<class RandomGenerator>
    void KSSolver<ScalarType, XCProvider>::initWaveFunc(RandomGenerator& gen) {
        {
            auto& waveUp = orbits[SpinState::Up][0];
            waveUp.flatten().random_normal(gen);
        }
        if constexpr (isSpinPolarized) {
            auto& waveDown = orbits[SpinState::Down][0];
            waveDown.flatten().random_normal(gen);
        }
    }

    template<class ScalarType, class XCProvider>
    template<class RandomGenerator>
    void KSSolver<ScalarType, XCProvider>::initWaveFunc(const SpinPair<KSOrbit, isSpinPolarized>& initial, RandomGenerator& gen) {
        const ssize_t dimX = std::min(initial[SpinState::Up].getDimX(), orbits[SpinState::Up][0].getDimX());
        const ssize_t dimY = std::min(initial[SpinState::Up].getDimY(), orbits[SpinState::Up][0].getDimY());
        const ssize_t dimZ = std::min(initial[SpinState::Up].getDimZ(), orbits[SpinState::Up][0].getDimZ());
        {
            auto& waveUp = orbits[SpinState::Up][0];
            waveUp.flatten().random_normal(gen);
            for (ssize_t x = -dimX; x <= dimX; ++x)
                for (ssize_t y = -dimY; y <= dimY; ++y)
                    for (ssize_t z = 0; z <= dimZ; ++z)
                        waveUp(x, y, z) = initial[SpinState::Up](x, y, z);
        }
        if constexpr (isSpinPolarized) {
            auto& waveDown = orbits[SpinState::Down][0];
            waveDown.flatten().random_normal(gen);
            for (ssize_t x = -dimX; x <= dimX; ++x)
                for (ssize_t y = -dimY; y <= dimY; ++y)
                    for (ssize_t z = 0; z <= dimZ; ++z)
                        waveDown(x, y, z) = initial[SpinState::Down](x, y, z);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initExternalPot() {
        using GridType = KSpaceGrid<ComplexType>;
        using Index3D = typename GridType::Index3D;
        using VectorType = Vector<ScalarType, 3>;
        externalPot = GridType::makeGrid(cutEnergy, repCell.getLattice());
        externalPot.flatten() = ScalarType(0);
        const auto factorGrids = cell.makeStructureFactor(cutEnergy);

        size_t i = 0;
        const std::unordered_set<uint16_t> species = cell.getSpecies();
        for (uint16_t element : species) {
            const ScalarType coeff = ScalarType(-4 * M_PI * PhyConst<AU>::unitCharge * PhyConst<AU>::unitCharge) * getCharge(element) / cell.getVolume();
            const auto& grid = factorGrids[i];
            GridType::forReducedKIndexInGrid(externalPot, repCell.getLattice(), [this, coeff, &grid](VectorType k, Index3D index) {
                const ScalarType squaredNorm = k.squaredNorm();
                const bool isNotGammaPoint = squaredNorm > std::numeric_limits<ScalarType>::min();
                if (isNotGammaPoint)
                    externalPot(index) += coeff * grid(index) / squaredNorm;
            });
            i += 1;
        }
        externalPot(0, 0, 0) = ComplexType(0);
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::assembleH(Vector3D k) {
        using VectorType = typename KSpaceGrid<ComplexType>::VectorType;
        h[SpinState::Up] = ScalarType(0);
        if constexpr (isSpinPolarized)
            h[SpinState::Down] = ScalarType(0);

        size_t i = 0;
        KSpaceGrid<ComplexType>::forKInGrid(orbits[SpinState::Up][0], repCell.getLattice(), [this, &i, &k](VectorType K) {
            constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
            const ScalarType kinetic = (k + K).squaredNorm() * factor;
            h[SpinState::Up](i, i) = kinetic;
            if constexpr (isSpinPolarized)
                h[SpinState::Down](i, i) = kinetic;
            i += 1;
        });
        fillPotential();
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::fillPotential() {
        using SignedIndex3D = typename KSpaceGrid<ComplexType>::Index3D;
        using VectorType = typename KSpaceGrid<ComplexType>::VectorType;

        RSpaceGrid<ComplexType> kSpaceXC_up{};
        RSpaceGrid<ComplexType> kSpaceXC_down{};
        RSpaceGrid<ComplexType> kSpaceDencity{};
        /* To kSpace */ {
            xcProvider.fill(density, xcPot);
            {
                fft.getRSpace().flatten() = xcPot[SpinState::Up].flatten();
                fft.transform();
                kSpaceXC_up = fft.getKSpace();
            }
            if constexpr (isSpinPolarized) {
                fft.getRSpace().flatten() = xcPot[SpinState::Down].flatten();
                fft.transform();
                kSpaceXC_down = fft.getKSpace();
            }
            fft.getRSpace().flatten() = density[SpinState::Up].flatten();
            fft.transform();
            kSpaceDencity = fft.getKSpace();
        }

        size_t row = 0;
        KSpaceGrid<ComplexType>::forKIndexInGrid(orbits[SpinState::Up][0], repCell.getLattice(),
            [this, &row, &kSpaceDencity, &kSpaceXC_up, &kSpaceXC_down](VectorType k1, SignedIndex3D index) {
                const ScalarType coeff = 1 / PhyConst<AU>::vacuumDielectric;
                const ScalarType repVolume = reciprocal(cell.getVolume());
                const auto any_wave = orbits[SpinState::Up][0];
                const ssize_t x1 = index[0];
                const ssize_t y1 = index[1];
                const ssize_t z1 = index[2];
                size_t col = row + 1;
                for (ssize_t x2 = x1; x2 < externalPot.getDimX(); ++x2) {
                    ssize_t y2 = x2 == x1 ? y1 : -externalPot.getDimY();
                    for (; y2 < externalPot.getDimY(); ++y2) {
                        ssize_t z2 = (x1 == x2 && y1 == y2) ? z1 + 1 : -externalPot.getDimZ();
                        for (; z2 < externalPot.getDimZ(); ++z2) {
                            SignedIndex3D delta{x1 - x2, y1 - y2, z1 - z2};
                            KSpaceGrid<ComplexType>::normalizeIndex(delta, externalPot.getDim());
                            const ssize_t deltaX = delta[0];
                            const ssize_t deltaY = delta[1];
                            const ssize_t deltaZ = delta[2];
                            const VectorType k2 = any_wave.getWaveVector(x2, y2, z2);
                            const VectorType deltaK = k1 - k2;

                            const ComplexType external = externalPot.calc(deltaX, deltaY, deltaZ);
                            const ComplexType hartree = kSpaceDencity(deltaX, deltaY, deltaZ) * coeff / deltaK.squaredNorm();
                            {
                                const ComplexType xc_up = kSpaceXC_up(deltaX, deltaY, deltaZ);
                                h[SpinState::Up](row, col) += (xc_up + hartree) * repVolume + external;
                            }
                            if constexpr (isSpinPolarized) {
                                const ComplexType xc_down = kSpaceXC_down(deltaX, deltaY, deltaZ);
                                h[SpinState::Down](row, col) += (xc_down + hartree) * repVolume + external;
                            }
                            col += 1;
                        }
                    }
                }
                row += 1;
            });
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateOrbits() {
        {
            auto& eigSolverUp = eigSolver[SpinState::Up];
            auto& orbits_up = orbits[SpinState::Up];
            const size_t orbitCount = orbits_up.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_up[i] = eigSolverUp.getEigenvectors().col(i);
        }
        if constexpr (isSpinPolarized) {
            auto& eigSolverDown = eigSolver[SpinState::Down];
            auto& orbits_down = orbits[SpinState::Down];
            const size_t orbitCount = orbits_down.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_down[i] = eigSolverDown.getEigenvectors().col(i);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateDensity(Vector3D k) {
        using UnsignedIndex3D = typename RSpaceGrid<ScalarType>::Index3D;
        using VectorType = Vector<ScalarType, 3>;
        chargeMixer.fetchInputDensity(density);

        const ScalarType inv_volume = reciprocal(cell.getVolume());
        auto& density_up = density[SpinState::Up];
        if constexpr (isSpinPolarized) {
            auto& density_down = density[SpinState::Down];
            auto density_compute_kernel = [this, k, &density_up, &density_down](VectorType pos, UnsignedIndex3D index) {
                ScalarType rho_up = ScalarType(0);
                for (const auto& orbit : orbits[SpinState::Up])
                    rho_up += orbit(k, pos).squaredNorm();
                density_up(index) = rho_up;

                ScalarType rho_down = ScalarType(0);
                for (const auto& orbit : orbits[SpinState::Down])
                    rho_down += orbit(k, pos).squaredNorm();
                density_down(index) = rho_down;
            };
            RSpaceGrid<ScalarType>::template forPointIndexInGrid<true, decltype(density_compute_kernel)>(
                density_up, cell.getLattice(), density_compute_kernel);
            density_up.flatten() *= inv_volume;
            density_down.flatten() *= inv_volume;
            /* Change format */ {
                auto rho = density_up.flatten();
                auto zeta = density_down.flatten();
                rho += zeta;
                zeta = divide(rho - zeta * ScalarType(2), rho);
            }
        }
        else {
            auto density_compute_kernel = [this, k, &density_up](VectorType pos, UnsignedIndex3D index) {
                auto rho_up = ScalarType(0);
                const auto& orbitsUp = orbits[SpinState::Up];
                size_t i = 0;
                for (; i < orbitsUp.getLength() - 1; ++i)
                    rho_up += orbitsUp[i](k, pos).squaredNorm() * ScalarType(2);
                rho_up += orbitsUp[i](k, pos).squaredNorm() * ScalarType(int(cell.getElectronCount() % 2U == 0) + 1);
                density_up(index) = rho_up;
            };
            RSpaceGrid<ScalarType>::template forPointIndexInGrid<true, decltype(density_compute_kernel)>(
                density_up, cell.getLattice(), density_compute_kernel);
            density_up.flatten() *= inv_volume;
        }
        chargeMixer.updateOutputDensity(density);
    }

    template<class ScalarType, class XCProvider>
    typename RSpaceGrid<ScalarType>::Index3D
    KSSolver<ScalarType, XCProvider>::makeFFTDim(ScalarType fftDimFactor) {
        if (fftDimFactor < ScalarType(2))
            throw std::invalid_argument("[Error]: FFT Dimension is too small");
        const auto dim = KSpaceGrid<ComplexType>::makeGridDim(cutEnergy * square(fftDimFactor), repCell.getLattice());
        typename RSpaceGrid<ScalarType>::Index3D result{};
        result[0] = dim[0] * 2 + 1;
        result[1] = dim[1] * 2 + 1;
        result[2] = dim[2] * 2 + 1;
        return result;
    }
}
