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
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/ElectronicStructure/ReciprocalCell.h"
#include "Ewald.h"
#include "PWBaseWave.h"
#include "Grid3D.h"

namespace Physica::Core::Internal {
    /**
     * \class KSSolverImpl extract \tparam isSpinPolarized from \tparam XCProvider
     */
    template<class ScalarType, class XCProvider, bool isSpinPolarized> class KSSolverImpl;

    template<class ScalarType, class XCProvider>
    class KSSolverImpl<ScalarType, XCProvider, true> {
        static_assert(XCProvider::isSpinPolarized == true);
    protected:
        using ComplexType = ComplexScalar<ScalarType>;
        using BandType = BandGrid<ScalarType, true>;
        using Vector3D = Vector<ScalarType, 3>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using Hamilton = std::pair<HermiteMatrix, HermiteMatrix>;
        using KSOrbit = PWBaseWave<ScalarType>;
        using KSOrbitArray = Utils::Array<KSOrbit>;
        using KSOrbits = std::pair<KSOrbitArray, KSOrbitArray>;
        using MatrixType = DenseMatrix<ComplexType>;
        using EigenSolverType = std::pair<EigenSolver<MatrixType>, EigenSolver<MatrixType>>;
        using UncenteredGrid = Grid3D<ScalarType, false>;
        using UnsignedDim = typename UncenteredGrid::Dim;
        using DensityType = std::pair<UncenteredGrid, UncenteredGrid>;
        using PotType = std::pair<UncenteredGrid, UncenteredGrid>;
        using CenteredGrid = Grid3D<ComplexType, true>;
        using SignedDim = typename CenteredGrid::Dim;

        constexpr static size_t DIISBufferSize = 3;
        using DensityRecord = Utils::Array<DensityType, DIISBufferSize>;
        using DIISBuffer = Utils::Array<UncenteredGrid, DIISBufferSize - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, DIISBufferSize, DIISBufferSize>;

        CrystalCell cell;
        ReciprocalCell repCell;
        ScalarType cutEnergy;
        KSOrbits orbits;
        EigenSolverType eigSolver;
        Hamilton h;
        BandType band;
        DensityRecord densityRecord;
        PotType xcPot;
        FFT<ScalarType, 3>* fft_xc_up;
        FFT<ScalarType, 3>* fft_xc_down;
        FFT<ScalarType, 3>* fft_hartree;
        XCProvider xcProvider;
        size_t iteration;
    public:
        KSSolverImpl(CrystalCell cell_, ScalarType cutEnergy_, BandType band_, UnsignedDim potGridDim);
        KSSolverImpl(const KSSolverImpl&) = delete;
        KSSolverImpl(KSSolverImpl&&) noexcept = delete;
        ~KSSolverImpl();
        /* Operators */
        KSSolverImpl& operator=(const KSSolverImpl& base) = delete;
        KSSolverImpl& operator=(KSSolverImpl&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte, const CenteredGrid& externalPot);
    protected:
        /* Operations */
        void initDensity();
        void assembleH(Vector3D k, const CenteredGrid& externalPot);
        void fillPotential(const CenteredGrid& externalPot);
        void updateOrbits();
        void updateDensity();
        void preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat);
        void DIISExtrapolation(DIISMatrix& diisMat);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return orbits.first[0].getPlainWaveCount(); }
        [[nodiscard]] size_t getDimX() const noexcept { return xcPot.first.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return xcPot.first.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return xcPot.first.getDimZ(); }
        [[nodiscard]] auto getDim() const noexcept { return xcPot.first.getDim(); }
        [[nodiscard]] size_t getSize() const noexcept { return xcProvider.getBufferSize(); }
        [[nodiscard]] size_t numOrbitToSolve() const { return (cell.getElectronCount() + 1) / 2; }
        [[nodiscard]] DensityType& currentDensity() { return *densityRecord.rbegin(); }
        [[nodiscard]] auto dimToPos(UnsignedDim dim) const noexcept { return xcPot.first.dimToPos(dim); }
        [[nodiscard]] SignedDim indexToSignedDim(size_t index) const noexcept { return orbits.first[0].indexToDim(index); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(SignedDim dim) const noexcept { return orbits.first[0].getWaveVector(dim); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(size_t index) const noexcept { return orbits.first[0].getWaveVector(index); }
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType, class XCProvider>
    KSSolverImpl<ScalarType, XCProvider, true>::KSSolverImpl(CrystalCell cell_,
                                                             ScalarType cutEnergy_,
                                                             BandType band_,
                                                             UnsignedDim potGridDim)
            : cell(std::move(cell_))
            , repCell(cell_.reciprocal())
            , cutEnergy(std::move(cutEnergy_))
            , band(std::move(band_))
            , densityRecord(DIISBufferSize, std::make_pair(UncenteredGrid(cell_.getLattice(), potGridDim),
                                                           UncenteredGrid(cell_.getLattice(), potGridDim)))
            , xcPot(std::make_pair(UncenteredGrid(cell_.getLattice(), potGridDim),
                                   UncenteredGrid(cell_.getLattice(), potGridDim)))
            , xcProvider(std::get<0>(potGridDim) * std::get<1>(potGridDim) * std::get<2>(potGridDim))
            , iteration(0) {
        const size_t electronCount = cell.getElectronCount();
        orbits = std::make_pair(KSOrbitArray((electronCount + 1) / 2, KSOrbit(cutEnergy, repCell.getLattice())),
                                KSOrbitArray(electronCount / 2, KSOrbit(cutEnergy, repCell.getLattice())));

        const size_t plainWaveCount = getPlainWaveCount();
        h = std::make_pair(HermiteMatrix(plainWaveCount), HermiteMatrix(plainWaveCount));
        eigSolver = std::make_pair(EigenSolver<MatrixType>(plainWaveCount), EigenSolver<MatrixType>(plainWaveCount));

        const Utils::Array<size_t, 3> fftGrid{std::get<0>(potGridDim), std::get<1>(potGridDim), std::get<2>(potGridDim)};
        const Utils::Array<ScalarType, 3> fftDeltaTs{ScalarType(cell.getLattice().row(0).norm()) / ScalarType(fftGrid[0] - 1),
                                                     ScalarType(cell.getLattice().row(1).norm()) / ScalarType(fftGrid[1] - 1),
                                                     ScalarType(cell.getLattice().row(2).norm()) / ScalarType(fftGrid[2] - 1)};
        fft_xc_up = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);
        fft_xc_down = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);
        fft_hartree = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);

        initDensity();
    }

    template<class ScalarType, class XCProvider>
    KSSolverImpl<ScalarType, XCProvider, true>::~KSSolverImpl() {
        delete fft_xc_up;
        delete fft_xc_down;
        delete fft_hartree;
    }

    template<class ScalarType, class XCProvider>
    bool KSSolverImpl<ScalarType, XCProvider, true>::solve(const ScalarType& criteria, size_t maxIte, const CenteredGrid& externalPot) {
        auto densityResiduals = DIISBuffer(DIISBufferSize - 1, UncenteredGrid(cell.getLattice(), getDimX(), getDimY(), getDimZ()));
        auto diisMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -ScalarType::One());
        diisMat(0, 0) = ScalarType::Zero();

        auto& eigSolverUp = eigSolver.first;
        auto& eigSolverDown = eigSolver.second;

        for (auto& kPoint : band.getKPoints()) {
            iteration = 0;
            while (true) {
                assembleH(kPoint.getPos(), externalPot);
                eigSolverUp.compute(h.first, true);
                eigSolverDown.compute(h.second, true);
                eigSolverUp.sort();
                eigSolverDown.sort();

                updateOrbits();
                updateDensity();

                if (iteration != 0) {
                    const auto& delta_rho = (*densityResiduals.crbegin()).asVector();
                    const auto& rho = currentDensity().first.asVector();
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
            kPoint.setEigInfo(eigSolverUp, eigSolverDown);
        }
        return true;
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, true>::initDensity() {
        const ScalarType averageDensity = ScalarType(cell.getElectronCount()) / cell.getVolume();
        auto& pair1 = currentDensity();
        auto& rho1 = pair1.first.asVector();
        rho1 = averageDensity;
        auto& zeta1 = pair1.second.asVector();
        zeta1 = ScalarType::Zero();

        auto& pair2 = densityRecord[densityRecord.getLength() - 2];
        auto& rho2 = pair2.first.asVector();
        rho2 = ScalarType::Zero();
        auto& zeta2 = pair2.second.asVector();
        zeta2 = ScalarType::Zero();
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, true>::assembleH(Vector3D k, const CenteredGrid& externalPot) {
        auto& h_up = h.first;
        auto& h_down = h.second;
        h_up = ScalarType::Zero();
        h_down = ScalarType::Zero();
        /* fill kinetic */ {
            const size_t order = h_up.getRow();
            for (size_t i = 0; i < order; ++i) {
                const ScalarType temp = ScalarType((k + getWaveVector(i)).squaredNorm()) * ScalarType(0.5);
                h_up(i, i) += temp;
                h_down(i, i) += temp;
            }
        }
        fillPotential(externalPot);
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, true>::fillPotential(const CenteredGrid& externalPot) {
        using VectorType = Vector<ScalarType, 3>;
        xcProvider.fill(currentDensity(), xcPot);
        fft_xc_up->transform(xcPot.first.asVector());
        fft_xc_down->transform(xcPot.second.asVector());
        fft_hartree->transform(currentDensity().first.asVector());

        const ScalarType factor = reciprocal(ScalarType(2 * M_PI));
        const auto fft_nomalizer = reciprocal(ScalarType(getSize()));
        const ScalarType factor1 = ScalarType(4 * M_PI) / cell.getVolume() * fft_nomalizer;

        auto& h_up = h.first;
        auto& h_down = h.second;
        const size_t order = h_up.getRow();
        for (size_t i = 0; i < order; ++i) {
            const auto dim1 = indexToSignedDim(i);
            auto[x1, y1, z1] = dim1;
            const VectorType k1 = getWaveVector(dim1);
            for (size_t j = i; j < order; ++j) {
                const auto dim2 = indexToSignedDim(j);
                auto[x2, y2, z2] = dim2;
                const VectorType k2 = getWaveVector(dim2);
                const VectorType deltaK = k1 - k2;
                const VectorType k = deltaK * factor;

                const ComplexType xc_up = fft_xc_up->getFreqIntense(k) * fft_nomalizer;
                const ComplexType xc_down = fft_xc_down->getFreqIntense(k) * fft_nomalizer;
                ComplexType hartree;
                if (i == j)
                    hartree = ComplexType::Zero();
                else
                    hartree = fft_hartree->getFreqIntense(k) * factor1 / deltaK.squaredNorm();
                const ComplexType external = externalPot(x1 - x2, y1 - y2, z1 - z2);

                h_up(i, j) += xc_up + hartree + external;
                h_down(i, j) += xc_down + hartree + external;
            }
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, true>::updateOrbits() {
        {
            auto& eigSolverUp = eigSolver.first;
            auto& orbits_up = orbits.first;
            const size_t orbitCount = orbits_up.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_up[i] = eigSolverUp.getRawEigenvectors().col(i); //TODO: Eigenvectors have been nomalized, find out the reason
        }
        {
            auto& eigSolverDown = eigSolver.second;
            auto& orbits_down = orbits.second;
            const size_t orbitCount = orbits_down.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_down[i] = eigSolverDown.getRawEigenvectors().col(i);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, true>::updateDensity() {
        /* Get density */ {
            auto& density_up = densityRecord[0].first;
            auto& density_down = densityRecord[0].second;
            auto& orbits_up = orbits.first;
            auto& orbits_down = orbits.second;
            auto[dimX, dimY, dimZ] = getDim();
            for (size_t i = 0; i < dimX; ++i) {
                for (size_t j = 0; j < dimY; ++j) {
                    for (size_t k = 0; k < dimZ; ++k) {
                        const auto pos = dimToPos({i, j, k});
                        auto rho_up = ScalarType::Zero();
                        for (size_t index = 0; index < orbits_up.getLength(); ++index)
                            rho_up += orbits_up[index](pos).squaredNorm();
                        density_up(i, j, k) = rho_up;

                        auto rho_down = ScalarType::Zero();
                        for (size_t index = 0; index < orbits_down.getLength(); ++index)
                            rho_down += orbits_down[index](pos).squaredNorm();
                        density_down(i, j, k) = rho_down;
                    }
                }
            }
            const ScalarType inv_volume = reciprocal(cell.getVolume());
            density_up.asVector() *= inv_volume;
            density_down.asVector() *= inv_volume;
        }
        /* Change format */ {
            auto& rho = densityRecord[0].first.asVector();
            auto& zeta = densityRecord[0].second.asVector();
            rho += zeta;
            zeta = divide(rho - zeta * ScalarType::Two(), rho);
        }

        for (size_t i = 0; i < densityRecord.getLength() - 1; ++i)
            swap(densityRecord[i], densityRecord[i + 1]);
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, true>::preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat) {
        /* Update residuals */ {
            const auto& rho_new = currentDensity().first.asVector();
            const auto& rho_old = densityRecord[densityRecord.getLength() - 2].first.asVector();
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
    void KSSolverImpl<ScalarType, XCProvider, true>::DIISExtrapolation(DIISMatrix& diisMat) {
        Vector<ScalarType, DIISBufferSize> x{};
        /* Solve linear equation */ {
            Vector<ScalarType, DIISBufferSize> b = Vector<ScalarType, DIISBufferSize>(DIISBufferSize, ScalarType::Zero());
            b[0] = -ScalarType::One();
            const DIISMatrix inv_A = diisMat.inverse();
            x = inv_A * b;
        }

        auto& new_rho = currentDensity().first.asVector();
        new_rho = ScalarType::Zero();
        for (size_t i = 1; i < x.getLength(); ++i) {
            const auto& rho = densityRecord[i - 1].first.asVector();
            new_rho += rho * x[i];
        }
    }

    template<class ScalarType, class XCProvider>
    class KSSolverImpl<ScalarType, XCProvider, false> {
        static_assert(XCProvider::isSpinPolarized == false);
    protected:
        using ComplexType = ComplexScalar<ScalarType>;
        using BandType = BandGrid<ScalarType, false>;
        using Vector3D = Vector<ScalarType, 3>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using Hamilton = HermiteMatrix;
        using KSOrbit = PWBaseWave<ScalarType>;
        using KSOrbitArray = Utils::Array<KSOrbit>;
        using KSOrbits = KSOrbitArray;
        using MatrixType = DenseMatrix<ComplexType>;
        using EigenSolverType = EigenSolver<MatrixType>;
        using UncenteredGrid = Grid3D<ScalarType, false>;
        using UnsignedDim = typename UncenteredGrid::Dim;
        using DensityType = UncenteredGrid;
        using PotType = UncenteredGrid;
        using CenteredGrid = Grid3D<ComplexType, true>;
        using SignedDim = typename CenteredGrid::Dim;

        constexpr static size_t DIISBufferSize = 3;
        using DensityRecord = Utils::Array<DensityType, DIISBufferSize>;
        using DIISBuffer = Utils::Array<UncenteredGrid, DIISBufferSize - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, DIISBufferSize, DIISBufferSize>;

        CrystalCell cell;
        ReciprocalCell repCell;
        ScalarType cutEnergy;
        KSOrbits orbits;
        EigenSolverType eigSolver;
        Hamilton h;
        BandType band;
        DensityRecord densityRecord;
        PotType xcPot;
        FFT<ScalarType, 3>* fft_xc;
        FFT<ScalarType, 3>* fft_hartree;
        XCProvider xcProvider;
        size_t iteration;
    public:
        KSSolverImpl(CrystalCell cell_, ScalarType cutEnergy_, BandType band_, UnsignedDim potGridDim);
        KSSolverImpl(const KSSolverImpl&) = delete;
        KSSolverImpl(KSSolverImpl&&) noexcept = delete;
        ~KSSolverImpl();
        /* Operators */
        KSSolverImpl& operator=(const KSSolverImpl& base) = delete;
        KSSolverImpl& operator=(KSSolverImpl&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte, const CenteredGrid& externalPot);
    protected:
        /* Operations */
        void initDensity();
        void assembleH(Vector3D k, const CenteredGrid& externalPot);
        void fillPotential(const CenteredGrid& externalPot);
        void updateOrbits();
        void updateDensity();
        void preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat);
        void DIISExtrapolation(DIISMatrix& diisMat);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return orbits[0].getPlainWaveCount(); }
        [[nodiscard]] size_t getDimX() const noexcept { return xcPot.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return xcPot.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return xcPot.getDimZ(); }
        [[nodiscard]] auto getDim() const noexcept { return xcPot.getDim(); }
        [[nodiscard]] size_t getSize() const noexcept { return xcProvider.getBufferSize(); }
        [[nodiscard]] size_t numOrbitToSolve() const { return (cell.getElectronCount() + 1) / 2; }
        [[nodiscard]] DensityType& currentDensity() { return *densityRecord.rbegin(); }
        [[nodiscard]] auto dimToPos(UnsignedDim dim) const noexcept { return xcPot.dimToPos(dim); }
        [[nodiscard]] SignedDim indexToSignedDim(size_t index) const noexcept { return orbits[0].indexToDim(index); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(SignedDim dim) const noexcept { return orbits[0].getWaveVector(dim); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(size_t index) const noexcept { return orbits[0].getWaveVector(index); }
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType, class XCProvider>
    KSSolverImpl<ScalarType, XCProvider, false>::KSSolverImpl(CrystalCell cell_,
                                                              ScalarType cutEnergy_,
                                                              BandType band_,
                                                              UnsignedDim potGridDim)
            : cell(std::move(cell_))
            , repCell(cell_.reciprocal())
            , cutEnergy(std::move(cutEnergy_))
            , band(std::move(band_))
            , densityRecord(DIISBufferSize, UncenteredGrid(cell_.getLattice(), potGridDim))
            , xcPot(UncenteredGrid(cell_.getLattice(), potGridDim))
            , xcProvider(std::get<0>(potGridDim) * std::get<1>(potGridDim) * std::get<2>(potGridDim))
            , iteration(0) {
        const size_t electronCount = cell.getElectronCount();
        orbits = KSOrbitArray(electronCount / 2, KSOrbit(cutEnergy, repCell.getLattice()));

        const size_t plainWaveCount = getPlainWaveCount();
        h = HermiteMatrix(plainWaveCount);
        eigSolver = EigenSolver<MatrixType>(plainWaveCount);

        const Utils::Array<size_t, 3> fftGrid{std::get<0>(potGridDim), std::get<1>(potGridDim), std::get<2>(potGridDim)};
        const Utils::Array<ScalarType, 3> fftDeltaTs{ScalarType(cell.getLattice().row(0).norm()) / ScalarType(fftGrid[0] - 1),
                                                     ScalarType(cell.getLattice().row(1).norm()) / ScalarType(fftGrid[1] - 1),
                                                     ScalarType(cell.getLattice().row(2).norm()) / ScalarType(fftGrid[2] - 1)};
        fft_xc = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);
        fft_hartree = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);

        initDensity();
    }

    template<class ScalarType, class XCProvider>
    KSSolverImpl<ScalarType, XCProvider, false>::~KSSolverImpl() {
        delete fft_xc;
        delete fft_hartree;
    }

    template<class ScalarType, class XCProvider>
    bool KSSolverImpl<ScalarType, XCProvider, false>::solve(const ScalarType& criteria, size_t maxIte, const CenteredGrid& externalPot) {
        auto densityResiduals = DIISBuffer(DIISBufferSize - 1, UncenteredGrid(cell.getLattice(), getDimX(), getDimY(), getDimZ()));
        auto diisMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -ScalarType::One());
        diisMat(0, 0) = ScalarType::Zero();

        for (auto& kPoint : band.getKPoints()) {
            iteration = 0;
            while (true) {
                assembleH(kPoint.getPos(), externalPot);
                eigSolver.compute(h, true);
                eigSolver.sort();

                updateOrbits();
                updateDensity();

                if (iteration != 0) {
                    const auto& delta_rho = (*densityResiduals.crbegin()).asVector();
                    const auto& rho = currentDensity().asVector();
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
            kPoint.setEigInfo(eigSolver);
        }
        return true;
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, false>::initDensity() {
        const ScalarType averageDensity = ScalarType(cell.getElectronCount()) / cell.getVolume();
        auto& rho1 = currentDensity().asVector();
        rho1 = averageDensity;

        auto& rho2 = densityRecord[densityRecord.getLength() - 2].asVector();
        rho2 = ScalarType::Zero();
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, false>::assembleH(Vector3D k, const CenteredGrid& externalPot) {
        h = ScalarType::Zero();
        /* fill kinetic */ {
            const size_t order = h.getRow();
            for (size_t i = 0; i < order; ++i) {
                const ScalarType temp = ScalarType((k + getWaveVector(i)).squaredNorm()) * ScalarType(0.5);
                h(i, i) += temp;
            }
        }
        fillPotential(externalPot);
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, false>::fillPotential(const CenteredGrid& externalPot) {
        using VectorType = Vector<ScalarType, 3>;
        xcProvider.fill(currentDensity(), xcPot);
        fft_xc->transform(xcPot.asVector());
        fft_hartree->transform(currentDensity().asVector());

        const ScalarType factor = reciprocal(ScalarType(2 * M_PI));
        const auto fft_nomalizer = reciprocal(ScalarType(getSize()));
        const ScalarType factor1 = ScalarType(4 * M_PI) / cell.getVolume() * fft_nomalizer;

        const size_t order = h.getRow();
        for (size_t i = 0; i < order; ++i) {
            const auto dim1 = indexToSignedDim(i);
            auto[x1, y1, z1] = dim1;
            const VectorType k1 = getWaveVector(dim1);
            for (size_t j = i; j < order; ++j) {
                const auto dim2 = indexToSignedDim(j);
                auto[x2, y2, z2] = dim2;
                const VectorType k2 = getWaveVector(dim2);
                const VectorType deltaK = k1 - k2;
                const VectorType k = deltaK * factor;

                const ComplexType xc = fft_xc->getFreqIntense(k) * fft_nomalizer;
                ComplexType hartree;
                if (i == j)
                    hartree = ComplexType::Zero();
                else
                    hartree = fft_hartree->getFreqIntense(k) * factor1 / deltaK.squaredNorm();
                const ComplexType external = externalPot(x1 - x2, y1 - y2, z1 - z2);

                h(i, j) += xc + hartree + external;
            }
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, false>::updateOrbits() {
        const size_t orbitCount = orbits.getLength();
        for (size_t i = 0; i < orbitCount; ++i)
            orbits[i] = eigSolver.getRawEigenvectors().col(i); //TODO: Eigenvectors have been nomalized, find out the reason
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, false>::updateDensity() {
        /* Get density */ {
            auto& density = densityRecord[0];
            auto& orbits_up = orbits;
            auto[dimX, dimY, dimZ] = getDim();
            for (size_t i = 0; i < dimX; ++i) {
                for (size_t j = 0; j < dimY; ++j) {
                    for (size_t k = 0; k < dimZ; ++k) {
                        const auto pos = dimToPos({i, j, k});
                        auto rho = ScalarType::Zero();
                        for (size_t index = 0; index < orbits_up.getLength(); ++index)
                            rho += orbits_up[index](pos).squaredNorm();
                        density(i, j, k) = rho * ScalarType::Two();
                    }
                }
            }
            const ScalarType inv_volume = reciprocal(cell.getVolume());
            density.asVector() *= inv_volume;
        }

        for (size_t i = 0; i < densityRecord.getLength() - 1; ++i)
            swap(densityRecord[i], densityRecord[i + 1]);
    }

    template<class ScalarType, class XCProvider>
    void KSSolverImpl<ScalarType, XCProvider, false>::preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat) {
        /* Update residuals */ {
            const auto& rho_new = currentDensity().asVector();
            const auto& rho_old = densityRecord[densityRecord.getLength() - 2].asVector();
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
    void KSSolverImpl<ScalarType, XCProvider, false>::DIISExtrapolation(DIISMatrix& diisMat) {
        Vector<ScalarType, DIISBufferSize> x{};
        /* Solve linear equation */ {
            Vector<ScalarType, DIISBufferSize> b = Vector<ScalarType, DIISBufferSize>(DIISBufferSize, ScalarType::Zero());
            b[0] = -ScalarType::One();
            const DIISMatrix inv_A = diisMat.inverse();
            x = inv_A * b;
        }

        auto& new_rho = currentDensity().asVector();
        new_rho = ScalarType::Zero();
        for (size_t i = 1; i < x.getLength(); ++i) {
            const auto& rho = densityRecord[i - 1].asVector();
            new_rho += rho * x[i];
        }
    }
}