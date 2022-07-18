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
#include "Physica/Core/Physics/ElectronicStructure/ReciprocalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/ReciprocalCell.h"
#include "BandGrid.h"
#include "PWBaseWave.h"
#include "Ewald.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType, class XCProvider>
    class KSSolver {
    public:
        constexpr static bool isSpinPolarized = XCProvider::isSpinPolarized;
        constexpr static size_t NumSpin = isSpinPolarized ? 2 : 1;
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using KSOrbit = PWBaseWave<ScalarType>;
        using KSOrbitArray = Utils::Array<KSOrbit>;
        using MatrixType = DenseMatrix<ComplexType>;
        using UncenteredGrid = Grid3D<ScalarType, false>;
        using UnsignedDim = typename UncenteredGrid::Dim;
        using CenteredGrid = Grid3D<ComplexType, true>;
        using SignedDim = typename CenteredGrid::Dim;
        using BandType = BandGrid<ScalarType, isSpinPolarized>;
        using KSOrbits = SpinPair<KSOrbitArray, isSpinPolarized>;
        using Hamilton = SpinPair<HermiteMatrix, isSpinPolarized>;
        using EigenSolverType = SpinPair<EigenSolver<MatrixType>, isSpinPolarized>;
        using DensityType = Utils::Array<UncenteredGrid, NumSpin>;
        using PotType = SpinPair<UncenteredGrid, isSpinPolarized>;
        using FFTxc = SpinPair<FFT<ScalarType, 3>, isSpinPolarized>;

        constexpr static size_t DIISBufferSize = 3;
        using DIISBuffer = Utils::Array<UncenteredGrid, DIISBufferSize - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Element, DIISBufferSize, DIISBufferSize>;
        using DensityRecord = Utils::Array<DensityType, DIISBufferSize>;
    protected:
        CrystalCell cell;
        ReciprocalCell repCell;
        ScalarType cutEnergy;
        size_t iteration;

        KSOrbits orbits;
        Hamilton h;
        BandType band;

        EigenSolverType eigSolver;
        DensityRecord densityRecord;
        CenteredGrid externalPot;
        PotType xcPot;
        XCProvider xcProvider;
        FFTxc fft_xc;
        FFT<ScalarType, 3>* fft_hartree;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergy_, ScalarType gridDimFactor, BandType band_);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver();
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return orbits[SpinState::Up][0].getPlainWaveCount(); }
        [[nodiscard]] const BandType& getBand() const noexcept { return band; }
    protected:
        /* Operations */
        Utils::Array<CenteredGrid> makeStructureFactor(ScalarType factorCutoff);
        UnsignedDim makePotGridDim(ScalarType gridDimFactor);
        void initExternalPot();
        void initDensity();
        void assembleH(Vector3D k, const CenteredGrid& externalPot);
        void fillPotential(const CenteredGrid& externalPot);
        void updateOrbits();
        void updateDensity();
        void preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat);
        void DIISExtrapolation(DIISMatrix& diisMat);
        /* Getters */
        [[nodiscard]] size_t numOrbitToSolve() const { return (cell.getElectronCount() + 1) / 2; }
        [[nodiscard]] SignedDim indexToSignedDim(size_t index) const noexcept { return orbits[SpinState::Up][0].indexToDim(index); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(SignedDim dim) const noexcept { return orbits[SpinState::Up][0].getWaveVector(dim); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(size_t index) const noexcept { return orbits[SpinState::Up][0].getWaveVector(index); }
        [[nodiscard]] DensityType& currentDensity() { return *densityRecord.rbegin(); }
        [[nodiscard]] size_t getDimX() const noexcept { return xcPot[SpinState::Up].getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return xcPot[SpinState::Up].getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return xcPot[SpinState::Up].getDimZ(); }
        [[nodiscard]] auto getDim() const noexcept { return xcPot[SpinState::Up].getDim(); }
        [[nodiscard]] auto dimToPos(UnsignedDim dim) const noexcept { return xcPot[SpinState::Up].dimToPos(dim); }
        [[nodiscard]] size_t getDensitySize() const noexcept { return xcProvider.getBufferSize(); }
        /* Static members */
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::KSSolver(CrystalCell cell_, ScalarType cutEnergy_, ScalarType gridDimFactor, BandType band_)
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
            const size_t plainWaveCount = getPlainWaveCount();
            h = Hamilton(NumSpin, plainWaveCount);
            eigSolver = EigenSolverType(NumSpin, plainWaveCount);
        }
        /* Mesh grid related */ {
            auto potGridDim = makePotGridDim(gridDimFactor);
            densityRecord = DensityRecord(DIISBufferSize, NumSpin, cell.getLattice(), potGridDim);
            xcProvider = XCProvider(std::get<0>(potGridDim) * std::get<1>(potGridDim) * std::get<2>(potGridDim));
            xcPot = PotType(NumSpin, cell.getLattice(), potGridDim);
            /* Allocate fft */ {
                const Utils::Array<size_t, 3> fftGrid{std::get<0>(potGridDim), std::get<1>(potGridDim), std::get<2>(potGridDim)};
                const Utils::Array<ScalarType, 3> fftDeltaTs{ScalarType(cell.getLattice().row(0).norm()) / ScalarType(fftGrid[0] - 1),
                                                            ScalarType(cell.getLattice().row(1).norm()) / ScalarType(fftGrid[1] - 1),
                                                            ScalarType(cell.getLattice().row(2).norm()) / ScalarType(fftGrid[2] - 1)};
                fft_xc = FFTxc(NumSpin, fftGrid, fftDeltaTs);
                fft_hartree = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);
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
        auto densityResiduals = DIISBuffer(DIISBufferSize - 1, UncenteredGrid(cell.getLattice(), getDimX(), getDimY(), getDimZ()));
        auto diisMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -ScalarType::One());
        diisMat(0, 0) = ScalarType::Zero();

        for (auto& kPoint : band.getKPoints()) {
            iteration = 0;
            while (true) {
                assembleH(kPoint.getPos(), externalPot);
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
    Utils::Array<typename KSSolver<ScalarType, XCProvider>::CenteredGrid>
    KSSolver<ScalarType, XCProvider>::makeStructureFactor(ScalarType factorCutoff) {
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
                        const ScalarType phase = g * r.asVector();
                        temp += ComplexType(cos(phase), sin(phase));
                    }
                }
                factors.asVector()[i] = temp;
            }
            ++j;
        }
        return all_factors;
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initExternalPot() {
        const ScalarType factorCutoff = cutEnergy * 8;
        externalPot = CenteredGrid::gridFromCutEnergy(factorCutoff, repCell.getLattice());
        const Utils::Array<CenteredGrid> all_factors = makeStructureFactor(factorCutoff);
        const ScalarType factor1 = ScalarType(-4 * M_PI) / cell.getVolume();
        const std::unordered_set<uint16_t> species = cell.getSpecies();

        externalPot.asVector() = ScalarType::Zero();
        const size_t gridSize = externalPot.getSize();

        size_t j = 0;
        for (uint16_t element : species) {
            const CenteredGrid& factors = all_factors[j];
            for (size_t i = 0; i < gridSize; ++i)
                externalPot[i] += factor1 * getCharge(element) * factors[i] / factors.indexToPos(i).squaredNorm();
            ++j;
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
    void KSSolver<ScalarType, XCProvider>::assembleH(Vector3D k, const CenteredGrid& externalPot) {
        auto& h_up = h[SpinState::Up];
        h_up = ScalarType::Zero();
        if constexpr (isSpinPolarized)
            h[SpinState::Down] = ScalarType::Zero();
        /* fill kinetic */ {
            const size_t order = h_up.getRow();
            for (size_t i = 0; i < order; ++i) {
                const ScalarType temp = ScalarType((k + getWaveVector(i)).squaredNorm()) * ScalarType(0.5);
                h_up(i, i) += temp;
                if constexpr (isSpinPolarized)
                    h[SpinState::Down](i, i) += temp;
            }
        }
        fillPotential(externalPot);
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::fillPotential(const CenteredGrid& externalPot) {
        using VectorType = Vector<ScalarType, 3>;
        xcProvider.fill(currentDensity(), xcPot);
        fft_xc[SpinState::Up].transform(xcPot[SpinState::Up].asVector());
        if constexpr (isSpinPolarized)
            fft_xc[SpinState::Down].transform(xcPot[SpinState::Down].asVector());
        fft_hartree->transform(currentDensity()[0].asVector());

        const ScalarType factor = reciprocal(ScalarType(2 * M_PI));
        const auto fft_nomalizer = reciprocal(ScalarType(getDensitySize()));
        const ScalarType factor1 = ScalarType(4 * M_PI) / cell.getVolume() * fft_nomalizer;

        const size_t numPW = getPlainWaveCount();
        for (size_t i = 0; i < numPW; ++i) {
            const auto dim1 = indexToSignedDim(i);
            auto[x1, y1, z1] = dim1;
            const VectorType k1 = getWaveVector(dim1);
            for (size_t j = i; j < numPW; ++j) {
                const auto dim2 = indexToSignedDim(j);
                auto[x2, y2, z2] = dim2;
                const VectorType k2 = getWaveVector(dim2);
                const VectorType deltaK = k1 - k2;
                const VectorType freq = deltaK * factor;

                ComplexType hartree;
                if (i == j)
                    hartree = ComplexType::Zero();
                else
                    hartree = fft_hartree->getFreqIntense(freq) * factor1 / deltaK.squaredNorm();
                const ComplexType external = externalPot(x1 - x2, y1 - y2, z1 - z2);

                const ComplexType xc_up = fft_xc[SpinState::Up].getFreqIntense(freq) * fft_nomalizer;
                h[SpinState::Up](i, j) += xc_up + hartree + external;
                if constexpr (isSpinPolarized) {
                    const ComplexType xc_down = fft_xc[SpinState::Down].getFreqIntense(freq) * fft_nomalizer;
                    h[SpinState::Down](i, j) += xc_down + hartree + external;
                }
            }
        }
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
        const auto[dimX, dimY, dimZ] = getDim();
        const ScalarType inv_volume = reciprocal(cell.getVolume());

        const size_t numOccupiedUp = numOrbitToSolve();
        const auto& orbits_up = orbits[SpinState::Up];
        auto& density_up = densityRecord[0][0];
        ScalarType total_density_up = 0;

        if constexpr (isSpinPolarized) {
            const auto& orbits_down = orbits[SpinState::Down];
            const size_t numOccupiedDown = numOccupiedUp - (cell.getElectronCount() % 2 != 0);
            auto& density_down = densityRecord[0][1];
            ScalarType total_density_down = 0;

            for (size_t i = 0; i < dimX; ++i) {
                for (size_t j = 0; j < dimY; ++j) {
                    for (size_t k = 0; k < dimZ; ++k) {
                        const auto pos = dimToPos({i, j, k});
                        auto rho_up = ScalarType::Zero();
                        for (size_t index = 0; index < numOccupiedUp; ++index)
                            rho_up += orbits_up[index](pos).squaredNorm();
                        density_up(i, j, k) = rho_up;
                        total_density_up += rho_up;

                        auto rho_down = ScalarType::Zero();
                        for (size_t index = 0; index < numOccupiedDown; ++index)
                            rho_down += orbits_down[index](pos).squaredNorm();
                        density_down(i, j, k) = rho_down;
                        total_density_down += rho_down;
                    }
                }
            }
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
            for (size_t i = 0; i < dimX; ++i) {
                for (size_t j = 0; j < dimY; ++j) {
                    for (size_t k = 0; k < dimZ; ++k) {
                        const auto pos = dimToPos({i, j, k});
                        auto rho_up = ScalarType::Zero();
                        for (size_t index = 0; index < numOccupiedUp; ++index)
                            rho_up += orbits_up[index](pos).squaredNorm();
                        const ScalarType temp = rho_up * ScalarType::Two();
                        density_up(i, j, k) = temp;
                        total_density_up += temp;
                    }
                }
            }
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
    typename KSSolver<ScalarType, XCProvider>::UnsignedDim
    KSSolver<ScalarType, XCProvider>::makePotGridDim(ScalarType gridDimFactor) {
        assert(gridDimFactor >= ScalarType::Two()); //Nyquist theory requests
        const auto& lattice = cell.getLattice();
        const auto& repLattice = repCell.getLattice();
        const auto[pwDimX, pwDimY, pwDimZ] = UncenteredGrid::dimFromCutEnergy(cutEnergy, repLattice);
        const ScalarType factor = gridDimFactor / ScalarType(2 * M_PI);

        size_t potDimX, potDimY, potDimZ;
        potDimX = size_t((ScalarType(repLattice.row(0).norm() * lattice.row(0).norm()) * factor * pwDimX).getTrivial());
        potDimY = size_t((ScalarType(repLattice.row(1).norm() * lattice.row(1).norm()) * factor * pwDimY).getTrivial());
        potDimZ = size_t((ScalarType(repLattice.row(2).norm() * lattice.row(2).norm()) * factor * pwDimZ).getTrivial());
        return {potDimX, potDimY, potDimZ};
    }
}
