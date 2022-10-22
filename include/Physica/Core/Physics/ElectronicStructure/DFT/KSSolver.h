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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/EigenSolver.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Math/Transform/FFTGrid.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "BandGrid.h"
#include "PWBaseWave.h"
#include "Ewald.h"
#include "RSpaceGrid.h"

namespace Physica::Core {
    template<class ScalarType, class XCProvider>
    class KSSolver {
    public:
        constexpr static bool isSpinPolarized = XCProvider::isSpinPolarized;
        constexpr static size_t NumSpin = isSpinPolarized ? 2 : 1;
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using RepCellType = ReciprocalCell<typename CrystalCell::ScalarType>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using KSOrbit = PWBaseWave<ScalarType>;
        using KSOrbitArray = Utils::Array<KSOrbit>;
        using MatrixType = DenseMatrix<ComplexType>;
        using BandType = BandGrid<ScalarType, isSpinPolarized>;
        using KSOrbits = SpinPair<KSOrbitArray, isSpinPolarized>;
        using Hamilton = SpinPair<HermiteMatrix, isSpinPolarized>;
        using EigenSolverType = SpinPair<EigenSolver<MatrixType>, isSpinPolarized>;
        using DensityType = Utils::Array<RSpaceGrid<ScalarType>, NumSpin>;
        using PotType = SpinPair<RSpaceGrid<ScalarType>, isSpinPolarized>;
        using FFTxc = SpinPair<FFT<ScalarType, 3>, isSpinPolarized>;

        constexpr static size_t DIISBufferSize = 3;
        using DIISBuffer = Utils::Array<RSpaceGrid<ScalarType>, DIISBufferSize - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Element, DIISBufferSize, DIISBufferSize>;
        using DensityRecord = Utils::Array<DensityType, DIISBufferSize>;
    protected:
        CrystalCell cell;
        RepCellType repCell;
        ScalarType cutEnergy;
        size_t iteration;

        KSOrbits orbits;
        Hamilton h;
        BandType band;

        EigenSolverType eigSolver;
        DensityRecord densityRecord;
        KSpaceGrid<ComplexType> externalPot;
        PotType xcPot;
        XCProvider xcProvider;
        FFTxc fft_xc;
        FFT<ScalarType, 3>* fft_hartree;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergy_, ScalarType fftDimFactor, BandType band_);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver();
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] size_t getNumPlaneWave() const noexcept { return orbits[SpinState::Up][0].getNumPlaneWave(); }
        [[nodiscard]] const BandType& getBand() const noexcept { return band; }
        [[nodiscard]] typename RSpaceGrid<ScalarType>::Index3D getFFTDim() const noexcept { return xcPot[SpinState::Up].getDim(); }
    protected:
        /* Operations */
        typename RSpaceGrid<ScalarType>::Index3D makeFFTDim(ScalarType fftDimFactor);
        void initExternalPot();
        void initDensity();
        void assembleH(Vector3D k);
        void fillPotential();
        void updateOrbits();
        void updateDensity();
        void preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat);
        void DIISExtrapolation(DIISMatrix& diisMat);
        /* Getters */
        [[nodiscard]] size_t numOrbitToSolve() const { return (cell.getElectronCount() + 1) / 2; }
        [[nodiscard]] DensityType& currentDensity() { return *densityRecord.rbegin(); }
        [[nodiscard]] size_t getDensitySize() const noexcept { return xcProvider.getBufferSize(); }
        /* Static members */
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
        static inline void normalizeIndex(ssize_t& index, ssize_t range) noexcept;
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
            orbits[SpinState::Up] = KSOrbitArray((electronCount + 1) / 2, cutEnergy, repCell.getLattice());
            if constexpr (isSpinPolarized)
                orbits[SpinState::Down] = KSOrbitArray(electronCount / 2, cutEnergy, repCell.getLattice());
        }
        /* Matrix related */ {
            const size_t planeWaveCount = getNumPlaneWave();
            h = Hamilton(NumSpin, planeWaveCount);
            eigSolver = EigenSolverType(NumSpin, planeWaveCount);
        }
        /* Mesh grid related */ {
            const auto dim = makeFFTDim(fftDimFactor);
            densityRecord = DensityRecord(DIISBufferSize, NumSpin, dim[0], dim[1], dim[2]);
            xcProvider = XCProvider(dim[0] * dim[1] * dim[2]);
            xcPot = PotType(NumSpin, dim[0], dim[1], dim[2]);
            /* Allocate fft */ {
                Utils::Array<ScalarType, 3> fftDeltaTs{};
                for (int i = 0; i < 3; ++i)
                    fftDeltaTs[i] = ScalarType(2 * M_PI) / repCell.getLattice().row(i).norm();
                fft_xc = FFTxc(NumSpin, dim, fftDeltaTs);
                fft_hartree = new FFT<ScalarType, 3>(dim, fftDeltaTs);
            }
        }
        initDensity();
        initExternalPot();
    }

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::~KSSolver() {
        delete fft_hartree;
    }

    template<class ScalarType, class XCProvider>
    bool KSSolver<ScalarType, XCProvider>::solve(const ScalarType& criteria, size_t maxIte) {
        const auto fft_dim = getFFTDim();
        auto densityResiduals = DIISBuffer(DIISBufferSize - 1, fft_dim[0], fft_dim[1], fft_dim[2]);
        auto diisMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -ScalarType::One());
        diisMat(0, 0) = ScalarType::Zero();

        for (auto& kPoint : band.getKPointGrid()) {
            iteration = 0;
            while (true) {
                assembleH(kPoint.getPos());
                eigSolver[SpinState::Up].compute(h[SpinState::Up], true);
                eigSolver[SpinState::Up].sort();
                if constexpr (isSpinPolarized) {
                    eigSolver[SpinState::Down].compute(h[SpinState::Down], true);
                    eigSolver[SpinState::Down].sort();
                }
                updateOrbits();
                updateDensity();

                if (iteration != 0) {
                    const auto& delta_rho = (*densityResiduals.crbegin()).asVector();
                    const auto& rho = currentDensity()[0].asVector();
                    const ScalarType error = abs(divide(delta_rho, rho)).max();
                    const bool isConverged = error < criteria;
                    if (isConverged)
                        break;
                }

                preDIIS(densityResiduals, diisMat);
                const bool doDIIS = iteration != 0 && iteration % DIISBufferSize == 0;
                if (doDIIS)
                    DIISExtrapolation(diisMat);

                if (++iteration == maxIte)
                    throw BadConvergenceException();
            };
            kPoint.setBandEnergy(SpinState::Up, toRealVector(eigSolver[SpinState::Up].getEigenvalues()));
            if constexpr (isSpinPolarized)
                kPoint.setBandEnergy(SpinState::Down, toRealVector(eigSolver[SpinState::Down].getEigenvalues()));
        }
        return true;
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initExternalPot() {
        using GridType = KSpaceGrid<ComplexType>;
        using VectorType = Vector<ScalarType, 3>;
        externalPot = GridType::makeGrid(cutEnergy, repCell.getLattice());
        externalPot.asVector() = ScalarType::Zero();
        const auto factorGrids = cell.makeStructureFactor(cutEnergy);

        size_t i = 0;
        const std::unordered_set<uint16_t> species = cell.getSpecies();
        for (uint16_t element : species) {
            const ScalarType coeff = ScalarType(-4 * M_PI * PhyConst<AU>::unitCharge * PhyConst<AU>::unitCharge) * getCharge(element) / cell.getVolume();
            const auto& grid = factorGrids[i++];
            size_t j = 0;
            GridType::forReducedKInGrid(grid.getDim(), repCell.getLattice(), [this, coeff, &grid, &j](VectorType k) {
                const ScalarType squaredNorm = k.squaredNorm();
                const bool isNotGammaPoint = squaredNorm > std::numeric_limits<ScalarType>::min();
                if (isNotGammaPoint)
                    externalPot.asVector()[j] += coeff * grid.asVector()[j] / squaredNorm;
                j += 1;
            });
            ++i;
        }
        externalPot(0, 0, 0) = ComplexType::Zero();
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initDensity() {
        const ScalarType averageDensity = ScalarType(cell.getElectronCount()) / cell.getVolume();
        auto& density = currentDensity();
        auto& rho1 = density[0].asVector();
        rho1 = averageDensity;
        auto& last_density = densityRecord[densityRecord.getLength() - 2];
        auto& rho2 = last_density[0].asVector();
        rho2 = ScalarType::Zero();
        if constexpr (isSpinPolarized) {
            auto& zeta1 = density[1].asVector();
            zeta1 = ScalarType::Zero();
            auto& zeta2 = last_density[1].asVector();
            zeta2 = ScalarType::Zero();
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::assembleH(Vector3D k) {
        using VectorType = typename KSpaceGrid<ComplexType>::VectorType;
        h[SpinState::Up] = ScalarType::Zero();
        if constexpr (isSpinPolarized)
            h[SpinState::Down] = ScalarType::Zero();

        size_t i = 0;
        KSpaceGrid<ComplexType>::forKInGrid(orbits[SpinState::Up][0].getDim(), repCell.getLattice(), [this, &i, &k](VectorType K) {
            constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
            const ScalarType kinetic = (k + K).squaredNorm() * factor;
            h[SpinState::Up](i, i) += kinetic;
            if constexpr (isSpinPolarized)
                h[SpinState::Down](i, i) += kinetic;
            i += 1;
        });
        fillPotential();
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::fillPotential() {
        using SignedIndex3D = typename KSpaceGrid<ComplexType>::Index3D;
        using VectorType = typename KSpaceGrid<ComplexType>::VectorType;

        FFTGrid<ComplexType> kSpaceXC_up{};
        FFTGrid<ComplexType> kSpaceXC_down{};
        FFTGrid<ComplexType> kSpaceDencity{};
        /* To kSpace */ {
            xcProvider.fill(currentDensity(), xcPot);

            fft_xc[SpinState::Up].transform(xcPot[SpinState::Up].asVector());
            kSpaceXC_up = FFTGrid<ComplexType>(fft_xc[SpinState::Up]);
            if constexpr (isSpinPolarized) {
                fft_xc[SpinState::Down].transform(xcPot[SpinState::Down].asVector());
                kSpaceXC_down = FFTGrid<ComplexType>(fft_xc[SpinState::Up]);
            }

            fft_hartree->transform(currentDensity()[0].asVector());
            kSpaceDencity = FFTGrid<ComplexType>(*fft_hartree);
        }

        size_t row = 0;
        KSpaceGrid<ComplexType>::forKIndexInGrid(orbits[SpinState::Up][0].getDim(), repCell.getLattice(),
            [this, &row, &kSpaceDencity, &kSpaceXC_up, &kSpaceXC_down](VectorType k1, SignedIndex3D index) {
                const ScalarType coeff = ScalarType(4 * M_PI) / cell.getVolume();
                const auto any_wave = orbits[SpinState::Up][0];
                const ssize_t x1 = index[0];
                const ssize_t y1 = index[1];
                const ssize_t z1 = index[2];
                size_t col = row + 1;
                for (ssize_t x2 = x1; x2 < externalPot.getDimX(); ++x2) {
                    for (ssize_t y2 = y1; y2 < externalPot.getDimY(); ++y2) {
                        for (ssize_t z2 = z1 + 1; z2 < externalPot.getDimZ(); ++z2) {
                            ssize_t deltaX = x1 - x2;
                            ssize_t deltaY = y1 - y2;
                            ssize_t deltaZ = z1 - z2;
                            normalizeIndex(deltaX, externalPot.getDimX());
                            normalizeIndex(deltaY, externalPot.getDimY());
                            normalizeIndex(deltaZ, externalPot.getDimZ());
                            const VectorType k2 = any_wave.getWaveVector(x2, y2, z2);
                            const VectorType deltaK = k1 - k2;


                            const ComplexType external = externalPot(deltaX, deltaY, deltaZ);
                            const ComplexType hartree = kSpaceDencity(deltaX, deltaY, deltaZ) * coeff / deltaK.squaredNorm();

                            const ComplexType xc_up = kSpaceXC_up(deltaX, deltaY, deltaZ);
                            h[SpinState::Up](row, col) += xc_up + hartree + external;
                            if constexpr (isSpinPolarized) {
                                const ComplexType xc_down = kSpaceXC_down(deltaX, deltaY, deltaZ);
                                h[SpinState::Down](row, col) += xc_down + hartree + external;
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
                orbits_up[i] = eigSolverUp.getRawEigenvectors().col(i);
        }
        if constexpr (isSpinPolarized) {
            auto& eigSolverDown = eigSolver[SpinState::Down];
            auto& orbits_down = orbits[SpinState::Down];
            const size_t orbitCount = orbits_down.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_down[i] = eigSolverDown.getRawEigenvectors().col(i);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateDensity() {
        using UnsignedIndex3D = typename RSpaceGrid<ScalarType>::Index3D;
        using VectorType = Vector<ScalarType, 3>;
        const ScalarType inv_volume = reciprocal(cell.getVolume());

        const size_t numOccupiedUp = numOrbitToSolve();
        auto& density_up = densityRecord[0][0];
        ScalarType total_density_up = 0;

        if constexpr (isSpinPolarized) {
            const size_t numOccupiedDown = numOccupiedUp - (cell.getElectronCount() % 2 != 0);
            auto& density_down = densityRecord[0][1];
            ScalarType total_density_down = 0;

            RSpaceGrid<ScalarType>::forPointIndexInGrid(density_up, cell.getLattice(),
                [this, &density_up, &density_down, &total_density_up, &total_density_down, numOccupiedUp, numOccupiedDown](VectorType pos, UnsignedIndex3D index) {
                    ScalarType rho_up = ScalarType::Zero();
                    for (const auto& orbit : orbits[SpinState::Up])
                        rho_up += orbit(pos).squaredNorm();
                    density_up(index) = rho_up;
                    total_density_up += rho_up;

                    ScalarType rho_down = ScalarType::Zero();
                    for (const auto& orbit : orbits[SpinState::Down])
                        rho_down += orbit(pos).squaredNorm();
                    density_down(index) = rho_down;
                    total_density_down += rho_down;
                });
            density_up.asVector() *= ScalarType(numOccupiedUp * density_up.asVector().getLength()) / total_density_up * inv_volume;
            density_down.asVector() *= ScalarType(numOccupiedDown * density_up.asVector().getLength()) / total_density_down * inv_volume;
            /* Change format */ {
                auto& rho = density_up.asVector();
                auto& zeta = density_down.asVector();
                rho += zeta;
                zeta = divide(rho - zeta * ScalarType::Two(), rho);
            }
        }
        else {
            RSpaceGrid<ScalarType>::forPointIndexInGrid(density_up, cell.getLattice(),
                [this, &density_up, &total_density_up, numOccupiedUp](VectorType pos, UnsignedIndex3D index) {
                    auto rho_up = ScalarType::Zero();
                    for (const auto& orbit : orbits[SpinState::Up])
                        rho_up += orbit(pos).squaredNorm();
                    const ScalarType temp = rho_up * ScalarType::Two();
                    density_up(index) = temp;
                    total_density_up += temp;
                });
            density_up.asVector() *= ScalarType(2 * numOccupiedUp * density_up.asVector().getLength()) / total_density_up * inv_volume;
        }

        for (size_t i = 0; i < densityRecord.getLength() - 1; ++i)
            swap(densityRecord[i], densityRecord[i + 1]);
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat) {
        /* Update residuals */ {
            const auto& rho_new = currentDensity()[0].asVector();
            const auto& rho_old = densityRecord[densityRecord.getLength() - 2][0].asVector();
            residuals[0].asVector() = rho_new - rho_old;
            for (size_t i = 0; i < residuals.getLength() - 1; ++i)
                swap(residuals[i], residuals[i + 1]);
        }
        /* Construct equation */ {
            for (size_t i = 1; i < diisMat.getRow(); ++i) {
                for (size_t j = i; j < diisMat.getColumn(); ++j) {
                    ScalarType temp = residuals[i - 1].asVector() * residuals[j - 1].asVector();
                    diisMat(i, j) = temp;
                    diisMat(j, i) = temp;
                }
            }
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::DIISExtrapolation(DIISMatrix& diisMat) {
        Vector<ScalarType, DIISBufferSize> x{};
        /* Solve linear equation */ {
            auto b = Vector<ScalarType, DIISBufferSize>(DIISBufferSize, ScalarType::Zero());
            b[0] = -ScalarType::One();
            const DIISMatrix inv_A = diisMat.inverse();
            x = inv_A * b;
        }

        auto& new_rho = currentDensity()[0].asVector();
        new_rho = ScalarType::Zero();
        for (size_t i = 1; i < x.getLength(); ++i) {
            const auto& rho = densityRecord[i - 1][0].asVector();
            new_rho += rho * x[i];
        }
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

     template<class ScalarType, class XCProvider>
    inline void KSSolver<ScalarType, XCProvider>::normalizeIndex(ssize_t& index, ssize_t range) noexcept {
        if (index > range)
            index -= range;
        else if (index < -range)
            index += range;
    }
}
