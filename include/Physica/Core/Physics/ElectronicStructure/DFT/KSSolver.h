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
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/ReciprocalCell.h"
#include "Ewald.h"
#include "WaveFunction.h"
#include "BandGrid.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType, class XCProvider>
    class KSSolver {
        using ComplexType = ComplexScalar<ScalarType>;
        using KPoint = typename BandGrid<ScalarType>::KPoint;
        using Hamilton = DenseHermiteMatrix<ComplexType>;
        using HamiltonPair = std::pair<Hamilton, Hamilton>;
        using MatrixType = DenseMatrix<ComplexType>;
        using KSOrbit = WaveFunction<ScalarType>;
        using KSOrbits = Utils::Array<KSOrbit>;
        using KSOrbitPair = std::pair<KSOrbits, KSOrbits>;
        using UncenteredGrid = Grid3D<ScalarType, false>;
        using UnsignedDim = typename UncenteredGrid::Dim;
        using DensityPair = std::pair<UncenteredGrid, UncenteredGrid>;
        using PotPair = std::pair<UncenteredGrid, UncenteredGrid>;
        using CenteredGrid = Grid3D<ComplexType, true>;
        using SignedDim = typename CenteredGrid::Dim;

        constexpr static size_t DIISBufferSize = 3;
        using DensityRecord = Utils::Array<DensityPair, DIISBufferSize>;
        using DIISBuffer = Utils::Array<UncenteredGrid, DIISBufferSize - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, DIISBufferSize, DIISBufferSize>;

        CrystalCell cell;
        ReciprocalCell repCell;
        ScalarType cutEnergy;
        KSOrbitPair orbits;
        BandGrid<ScalarType> kPoints;
        DensityRecord densityRecord;
        PotPair xcPot;
        CenteredGrid externalPotGrid;
        FFT<ScalarType, 3>* fft_xc_up;
        FFT<ScalarType, 3>* fft_xc_down;
        FFT<ScalarType, 3>* fft_hartree;
        XCProvider xcProvider;
        size_t iteration;
    public:
        KSSolver(CrystalCell cell_, ScalarType cutEnergy_, BandGrid<ScalarType> kPoints_, size_t gridDimX_, size_t gridDimY_, size_t gridDimZ_);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver();
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] const BandGrid<ScalarType>& getBandGrid() const noexcept { return kPoints; }
    private:
        /* Operations */
        void initialize(size_t gridDimX, size_t gridDimY, size_t gridDimZ);
        void initDensity();
        void initExternalPot();
        void assembleH(KPoint k, HamiltonPair& hPair);
        void fillPotential(HamiltonPair& hPair);
        void updateOrbits(const EigenSolver<MatrixType>& eigenSolver_up, const EigenSolver<MatrixType>& eigenSolver_down);
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
        [[nodiscard]] DensityPair& currentDensity() { return *densityRecord.rbegin(); }
        [[nodiscard]] auto dimToPos(UnsignedDim dim) const noexcept { return xcPot.first.dimToPos(dim); }
        [[nodiscard]] SignedDim indexToSignedDim(size_t index) const noexcept { return orbits.first[0].indexToDim(index); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(SignedDim dim) const noexcept { return orbits.first[0].getWaveVector(dim); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(size_t index) const noexcept { return orbits.first[0].getWaveVector(index); }
        [[nodiscard]] Utils::Array<CenteredGrid> getStructureFactor(ScalarType factorCutoff);
        [[nodiscard]] static int16_t getCharge(uint16_t atomicNum) { return atomicNum; }
    };

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::KSSolver(CrystalCell cell_,
                                               ScalarType cutEnergy_,
                                               BandGrid<ScalarType> kPoints_,
                                               size_t gridDimX,
                                               size_t gridDimY,
                                               size_t gridDimZ)
            : cell(std::move(cell_))
            , repCell(cell_.reciprocal())
            , cutEnergy(std::move(cutEnergy_))
            , kPoints(std::move(kPoints_))
            , densityRecord(DIISBufferSize, std::make_pair(UncenteredGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ),
                                                           UncenteredGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ)))
            , xcPot(std::make_pair(UncenteredGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ),
                                   UncenteredGrid(cell_.getLattice(), gridDimX, gridDimY, gridDimZ)))
            , xcProvider(gridDimX * gridDimY * gridDimZ)
            , iteration(0) {
        initialize(gridDimX, gridDimY, gridDimZ);
    }

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::~KSSolver() {
        delete fft_xc_up;
        delete fft_xc_down;
        delete fft_hartree;
    }

    template<class ScalarType, class XCProvider>
    bool KSSolver<ScalarType, XCProvider>::solve(const ScalarType& criteria, size_t maxIte) {
        const size_t plainWaveCount = getPlainWaveCount();
        auto hamilton = std::make_pair(Hamilton(plainWaveCount), Hamilton(plainWaveCount));

        auto densityResiduals = DIISBuffer(DIISBufferSize - 1, UncenteredGrid(cell.getLattice(), getDimX(), getDimY(), getDimZ()));
        auto diisMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -ScalarType::One());
        diisMat(0, 0) = ScalarType::Zero();

        auto eigenSolver_up = EigenSolver<MatrixType>(plainWaveCount);
        auto eigenSolver_down = EigenSolver<MatrixType>(plainWaveCount);

        const size_t kPointCount = kPoints.getSize();
        for (size_t i = 0; i < kPointCount; ++i) {
            const KPoint kPoint = kPoints.indexToPos(i);
            iteration = 0;
            while (true) {
                assembleH(kPoint, hamilton);
                eigenSolver_up.compute(hamilton.first, true);
                eigenSolver_down.compute(hamilton.second, true);
                eigenSolver_up.sort();
                eigenSolver_down.sort();

                updateOrbits(eigenSolver_up, eigenSolver_down);
                updateDensity();

                if (iteration != 0) {
                    const auto& last_rho = (*densityResiduals.crbegin()).asVector();
                    const auto& rho = currentDensity().first.asVector();
                    const ScalarType delta_rho = abs(divide(last_rho, rho)).max();
                    const bool isConverged = delta_rho < criteria;
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
            kPoints[i].setRawData(eigenSolver_up, eigenSolver_down);
        }
        return true;
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initialize(size_t gridDimX, size_t gridDimY, size_t gridDimZ) {
        const size_t electronCount = cell.getElectronCount();
        orbits = std::make_pair(KSOrbits((electronCount + 1) / 2, KSOrbit(cutEnergy, repCell.getLattice())),
                                KSOrbits(electronCount / 2, KSOrbit(cutEnergy, repCell.getLattice())));

        const Utils::Array<size_t, 3> fftGrid{gridDimX, gridDimY, gridDimZ};
        const Utils::Array<ScalarType, 3> fftDeltaTs{ScalarType(cell.getLattice().row(0).norm()) / ScalarType(gridDimX - 1),
                                                     ScalarType(cell.getLattice().row(1).norm()) / ScalarType(gridDimY - 1),
                                                     ScalarType(cell.getLattice().row(2).norm()) / ScalarType(gridDimZ - 1)};
        fft_xc_up = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);
        fft_xc_down = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);
        fft_hartree = new FFT<ScalarType, 3>(fftGrid, fftDeltaTs);

        initExternalPot();
        initDensity();
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initDensity() {
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
    void KSSolver<ScalarType, XCProvider>::initExternalPot() {
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

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::assembleH(KPoint k, HamiltonPair& hPair) {
        auto& h_up = hPair.first;
        auto& h_down = hPair.second;
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
        fillPotential(hPair);
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::fillPotential(HamiltonPair& hPair) {
        using VectorType = Vector<ScalarType, 3>;
        xcProvider.fill(currentDensity(), xcPot);
        fft_xc_up->transform(xcPot.first.asVector());
        fft_xc_down->transform(xcPot.second.asVector());
        fft_hartree->transform(currentDensity().first.asVector());

        const ScalarType factor = reciprocal(ScalarType(2 * M_PI));
        const auto fft_nomalizer = reciprocal(ScalarType(getSize()));
        const ScalarType factor1 = ScalarType(4 * M_PI) / cell.getVolume() * fft_nomalizer;
        
        auto& h_up = hPair.first;
        auto& h_down = hPair.second;
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
                const ComplexType xc_down = fft_xc_up->getFreqIntense(k) * fft_nomalizer;
                ComplexType hartree;
                if (i == j)
                    hartree = ComplexType::Zero();
                else
                    hartree = fft_hartree->getFreqIntense(k) * factor1 / deltaK.squaredNorm();
                const ComplexType external = externalPotGrid(x1 - x2, y1 - y2, z1 - z2);

                h_up(i, j) += xc_up + hartree + external;
                h_down(i, j) += xc_down + hartree + external;
            }
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateOrbits(const EigenSolver<MatrixType>& eigenSolver_up, const EigenSolver<MatrixType>& eigenSolver_down) {
        {
            auto& orbits_up = orbits.first;
            const size_t orbitCount = orbits_up.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_up[i] = eigenSolver_up.getRawEigenvectors().col(i); //TODO: Eigenvectors have been nomalized, find out the reason
        }
        {
            auto& orbits_down = orbits.second;
            const size_t orbitCount = orbits_down.getLength();
            for (size_t i = 0; i < orbitCount; ++i)
                orbits_down[i] = eigenSolver_down.getRawEigenvectors().col(i);
        }
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::updateDensity() {
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
            std::swap(densityRecord[i], densityRecord[i + 1]);
    }

    template<class ScalarType, class XCProvider>
    typename Utils::Array<typename KSSolver<ScalarType, XCProvider>::CenteredGrid> KSSolver<ScalarType, XCProvider>::getStructureFactor(ScalarType factorCutoff) {
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

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::preDIIS(DIISBuffer& residuals, DIISMatrix& diisMat) {
        /* Update residuals */ {
            const auto& rho_new = currentDensity().first.asVector();
            const auto& rho_old = densityRecord[densityRecord.getLength() - 2].first.asVector();
            residuals[0].asVector() = rho_new - rho_old;
            for (size_t i = 0; i < residuals.getLength() - 1; ++i)
                std::swap(residuals[i], residuals[i + 1]);
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
}
