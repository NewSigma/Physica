/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "DensityGrid.h"

namespace Physica::Core {
    template<class ScalarType, bool IsSpinPolarized>
    class KSHamilton {
        static_assert(!IsSpinPolarized, "[Error]: Not implemented");
    public:
        using ComplexType = ComplexScalar<ScalarType>;
        using GridType = RSpaceGrid<ComplexType>;
        using StrucFactorArray = Utils::Array<GridType>;
        using LatticeMatrix = typename PeriodicCell<ScalarType, 3>::LatticeMatrix;
        using BasisType = PlainWaveBasis<ScalarType>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using Vector3D = Vector<ScalarType, 3>;
        using Index3D = typename GridBase::Index3D;
        using FFT3D = FFT<ScalarType, 3>;
    private:
        HermiteMatrix hamiltonH;
        CrystalCell<ScalarType> cell;
        LatticeMatrix repLatt;
        ScalarType cutEnergyPsi;
        ScalarType cutEnergyRho;
        Index3D basisDim;

        GridType kSpaceIonCoulomb;
        FFT3D fft_rho;
        FFT3D fft_xc;
    public:
        KSHamilton(CrystalCell<ScalarType> cell_, ScalarType cutEnergyPsi_, ScalarType cutEnergyRho_);
        KSHamilton(const KSHamilton&) = default;
        KSHamilton(KSHamilton&&) noexcept = default;
        ~KSHamilton() = default;
        /* Operators */
        KSHamilton& operator=(KSHamilton obj) noexcept { swap(obj); return *this; }
        const HermiteMatrix& operator[]([[maybe_unused]] SpinState spin) const noexcept { return hamiltonH; }
        /* Operations */
        void makeFreeHamilton(Vector3D kPoint);
        void makeNearFreeHamilton(Vector3D kPoint);
        void makeHamiltonWithoutXC(const DensityGrid<ScalarType, IsSpinPolarized>& densityRho, Vector3D kPoint);
        void makeHamiltonWithXC(
            const RSpaceGrid<ScalarType>& xcPot,
            const DensityGrid<ScalarType, IsSpinPolarized>& densityRho,
            Vector3D kPoint);

        void swap(KSHamilton& obj);
        /* Getters */
        [[nodiscard]] const CrystalCell<ScalarType>& getCrystalCell() const noexcept { return cell; }
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const LatticeMatrix& getRepLattice() const noexcept { return repLatt; }
        [[nodiscard]] ScalarType getCutEnergyPsi() const noexcept { return cutEnergyPsi; }
        [[nodiscard]] ScalarType getCutEnergyRho() const noexcept { return cutEnergyRho; }
        [[nodiscard]] size_t getNumBasis() const noexcept { return basisDim[0] * basisDim[1] * basisDim[2]; }
        [[nodiscard]] Index3D getFFTRSpaceSize() const noexcept { return fft_rho.getRSpaceSize(); }
    private:
        Index3D calcDeltaIndex(Index3D index1, Index3D index2, Index3D deltaDim) const;
        StrucFactorArray makeStructureFactor() const;
        void makeKSpaceIonCoulomb();
    };

    template<class ScalarType, bool IsSpinPolarized>
    KSHamilton<ScalarType, IsSpinPolarized>::KSHamilton(CrystalCell<ScalarType> cell_, ScalarType cutEnergyPsi_, ScalarType cutEnergyRho_)
            : cell(std::move(cell_))
            , cutEnergyPsi(cutEnergyPsi_)
            , cutEnergyRho(cutEnergyRho_) {
        assert(cutEnergyPsi * 4.0 < cutEnergyRho && "[Error]: Cut energy for charge density is too small");
        repLatt = cell.reciprocal().getLattice();
        basisDim = BasisType::makeGridDim(cutEnergyPsi, repLatt);
        const Index3D fftDim = BasisType::makeGridDim(cutEnergyRho, repLatt);
        fft_rho = FFT3D(fftDim, {1, 1, 1}, PlanFlag::Estimate);
        fft_xc = FFT3D(fftDim, {1, 1, 1}, PlanFlag::Estimate);
        makeKSpaceIonCoulomb();
    }

    template<class ScalarType, bool IsSpinPolarized>
    void KSHamilton<ScalarType, IsSpinPolarized>::swap(KSHamilton& __restrict obj) {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        hamiltonH.swap(obj.hamiltonH);
        cell.swap(obj.cell);
        repLatt.swap(obj.repLatt);
        cutEnergyPsi.swap(obj.cutEnergyPsi);
        cutEnergyRho.swap(obj.cutEnergyRho);
        basisDim.swap(obj.basisDim);
        fft_rho.swap(obj.fft_rho);
        fft_xc.swap(obj.fft_xc);
    }

    template<class ScalarType, bool IsSpinPolarized>
    void KSHamilton<ScalarType, IsSpinPolarized>::makeFreeHamilton(Vector3D kPoint) {
        hamiltonH = ScalarType(0);
        size_t i = 0;
        BasisType::forKInBasis(repLatt, basisDim, [this, kPoint, &i](Vector3D waveK) {
            constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
            const ScalarType kineticE = (kPoint + waveK).squaredNorm() * factor;
            hamiltonH(i, i) = kineticE;
            i += 1;
        });
    }

    template<class ScalarType, bool IsSpinPolarized>
    void KSHamilton<ScalarType, IsSpinPolarized>::makeNearFreeHamilton(Vector3D kPoint) {
        makeFreeHamilton(kPoint);

        GridBase::forIndexInGrid(basisDim, [this](Index3D index1) {
            const size_t r = PeriodIndex3D(basisDim, index1).toIndex1D();
            GridBase::forIndexInGrid(basisDim, [this, r, index1](Index3D index2) {
                const size_t c = PeriodIndex3D(basisDim, index2).toIndex1D();
                const Index3D delta = calcDeltaIndex(index1, index2, kSpaceIonCoulomb.getDim());
                hamiltonH(r, c) += kSpaceIonCoulomb(delta);
            });
        });
    }

    template<class ScalarType, bool IsSpinPolarized>
    void KSHamilton<ScalarType, IsSpinPolarized>::makeHamiltonWithoutXC(
            const DensityGrid<ScalarType, IsSpinPolarized>& densityRho,
            Vector3D kPoint) {
        makeNearFreeHamilton(kPoint);

        fft_rho.getRSpace().flatten() = densityRho.getTotalDensity().flatten();
        fft_rho.transform();
        auto& kSpaceDensity = fft_rho.getKSpace();
        const ScalarType repVolume = reciprocal(cell.getVolume());
        GridBase::forIndexInGrid(basisDim, [this, repVolume, &kSpaceDensity](Index3D index1) {
            const size_t r = PeriodIndex3D(basisDim, index1).toIndex1D();
            GridBase::forIndexInGrid(basisDim, [this, repVolume, r, index1, &kSpaceDensity](Index3D index2) {
                constexpr double factor = 1 / PhyConst<AU>::vacuumDielectric;
                const size_t c = PeriodIndex3D(basisDim, index2).toIndex1D();
                const Index3D delta = calcDeltaIndex(index1, index2, kSpaceDensity.getDim());
                const Vector3D waveG = BasisType::makeWaveVector(repLatt, delta, basisDim);
                hamiltonH(r, c) += ScalarType(factor) * repVolume * kSpaceDensity(delta) / waveG.squaredNorm();
            });
        });
    }

    template<class ScalarType, bool IsSpinPolarized>
    void KSHamilton<ScalarType, IsSpinPolarized>::makeHamiltonWithXC(
            const RSpaceGrid<ScalarType>& xcPot,
            const DensityGrid<ScalarType, IsSpinPolarized>& densityRho,
            Vector3D kPoint) {
        makeHamiltonWithoutXC(densityRho, kPoint);

        fft_xc.getRSpace().flatten() = xcPot.flatten();
        fft_xc.transform();
        auto& kSpaceXC = fft_xc.getKSpace();
        const ScalarType repVolume = reciprocal(cell.getVolume());
        GridBase::forIndexInGrid(basisDim, [this, repVolume, &kSpaceXC](Index3D index1) {
            const size_t r = PeriodIndex3D(basisDim, index1).toIndex1D();
            GridBase::forIndexInGrid(basisDim, [this, repVolume, r, index1, &kSpaceXC](Index3D index2) {
                const size_t c = PeriodIndex3D(basisDim, index2).toIndex1D();
                const Index3D delta = calcDeltaIndex(index1, index2, kSpaceXC.getDim());
                hamiltonH(r, c) += kSpaceXC(delta) * repVolume;
            });
        });
    }

    template<class ScalarType, bool IsSpinPolarized>
    typename KSHamilton<ScalarType, IsSpinPolarized>::Index3D
    KSHamilton<ScalarType, IsSpinPolarized>::calcDeltaIndex(Index3D index1, Index3D index2, Index3D deltaDim) const {
        for (int i = 0; i < 3; ++i) {
            assert(index1[i] < basisDim[i]);
            assert(index2[i] < basisDim[i]);
        }

        Index3D delta{};
        for (int i = 0; i < 3; ++i) {
            const ssize_t sIndex1 = index1[i] > basisDim[i] / 2 ? -ssize_t(basisDim[i] - index1[i]) : index1[i];
            const ssize_t sIndex2 = index2[i] > basisDim[i] / 2 ? -ssize_t(basisDim[i] - index2[i]) : index2[i];
            const ssize_t sIndex = sIndex1 - sIndex2;
            delta[i] = sIndex > 0 ? sIndex : sIndex + ssize_t(deltaDim[i]);
            assert(delta[i] >= 0);
            assert(delta[i] < deltaDim[i]);
        }
        return delta;
    }

    template<class ScalarType, bool IsSpinPolarized>
    typename KSHamilton<ScalarType, IsSpinPolarized>::StrucFactorArray
    KSHamilton<ScalarType, IsSpinPolarized>::makeStructureFactor() const {
        const auto species = cell.getSpecies();
        auto result = Utils::Array<GridType>(species.size(), kSpaceIonCoulomb.getDim());
        size_t i = 0;
        for (uint16_t element : species) {
            GridType& grid = result[i];
            i += 1;
            size_t j = 0;
            BasisType::forKInBasis(repLatt, grid.getDim(), [this, element, &j, &grid](Vector3D waveK) {
                auto factor = ComplexType(0);
                for (size_t ion = 0; ion < cell.getNumParticle(); ++ion) {
                    if (cell.getAtomicNumber(ion) == element) { //Optimize: We can use searching table method
                        const ScalarType phase = waveK * cell.getPos().row(ion).asVector();
                        factor += ComplexType::fromPhase(phase);
                    }
                }
                grid.flatten()[j] = factor;
                j += 1;
            });
        }
        return result;
    }

    template<class ScalarType, bool IsSpinPolarized>
    void KSHamilton<ScalarType, IsSpinPolarized>::makeKSpaceIonCoulomb() {
        {
            Index3D dim{};
            for (int i = 0; i < 3; ++i)
                dim[i] = basisDim[i] * 2 - 1;
            kSpaceIonCoulomb.resize(dim);
        }
        const std::unordered_set<uint16_t> species = cell.getSpecies();
        const auto strucFactorGrids = makeStructureFactor();
        const ScalarType factor1 = ScalarType(-4 * M_PI * PhyConst<AU>::unitCharge * PhyConst<AU>::unitCharge) / cell.getVolume();
        size_t i = 0;
        for (uint16_t element : species) {
            const ScalarType chargeZ = ScalarType(element);
            const ScalarType factor2 = factor1 * chargeZ;
            const auto& grid = strucFactorGrids[i];
            BasisType::forKIndexInBasis(repLatt, kSpaceIonCoulomb.getDim(), [this, factor2, &grid](Vector3D waveK, Index3D index) {
                const ScalarType squaredNorm = waveK.squaredNorm();
                const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                if (isNotGammaPoint) [[likely]]
                    kSpaceIonCoulomb(index) += factor2 * grid(index) / squaredNorm;
            });
            i += 1;
        }
        kSpaceIonCoulomb(0, 0, 0) = ScalarType(0);
    }
}
