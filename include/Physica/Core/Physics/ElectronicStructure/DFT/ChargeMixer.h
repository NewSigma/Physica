/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Math/Transform/FFT.h"
#include "DensityGrid.h"

namespace Physica::Core {
    /**
     * References:
     * [1] G. Kresse and J. Furthmüller, Efficient iterative schemes for ab initio total-energy calculations using a plane-wave basis set(https://doi.org/10.1103/physrevb.54.11169)
     */
    template<class ScalarType, bool isSpinPolarized>
    class ChargeMixer {
        constexpr static size_t DIISBufferSize = 5;
        constexpr static double amix = 0.8; //Refer to [1]
        constexpr static double bmix = 0.8; //Refer to [1]
        constexpr static double amin = 0.4; //Refer to VASP wiki
        constexpr static double pulay_mix = 0.4;
        using LatticeMatrix = typename CrystalCell<ScalarType>::LatticeMatrix;
        using DensityType = DensityGrid<ScalarType, isSpinPolarized>;
        using DensityArray = Array<DensityType, DIISBufferSize>;
        using FFT3D = FFT<ScalarType, 3>;
        using Index3D = typename DensityType::Index3D;
    private:
        LatticeMatrix repLatt;
        DensityArray oldDensities;
        DensityArray residules;
        FFT3D fft;
        size_t mixIteration;
    public:
        ChargeMixer() = default;
        ChargeMixer(LatticeMatrix repLatt_, Index3D dim);
        ChargeMixer(const ChargeMixer&) = default;
        ChargeMixer(ChargeMixer&&) noexcept = default;
        ~ChargeMixer() = default;
        /* Operators */
        ChargeMixer& operator=(ChargeMixer mixer) noexcept;
        /* Operations */
        void fetchInputDensity(DensityType& input);
        void mix(const DensityType& output);
        void mix(size_t iteration, DensityType& result);
        void swap(ChargeMixer& __restrict mixer) noexcept;
        /* Getters */
        [[nodiscard]] const DensityType& getResidule() const noexcept { return residules[mixIteration]; }
    };

    template<class ScalarType, bool isSpinPolarized>
    ChargeMixer<ScalarType, isSpinPolarized>::ChargeMixer(LatticeMatrix repLatt_, Index3D dim)
            : repLatt(std::move(repLatt_))
            , oldDensities(DIISBufferSize, dim)
            , residules(DIISBufferSize, dim)
            , fft(dim, PlanFlag::Estimate)
            , mixIteration(0) {}

    template<class ScalarType, bool isSpinPolarized>
    ChargeMixer<ScalarType, isSpinPolarized>& ChargeMixer<ScalarType, isSpinPolarized>::operator=(ChargeMixer mixer) noexcept {
        swap(mixer);
        return *this;
    }

    template<class ScalarType, bool isSpinPolarized>
    void ChargeMixer<ScalarType, isSpinPolarized>::swap(ChargeMixer& __restrict mixer) noexcept {
        assert(this != &mixer && "[Error]: Self swap is likely a bug");
        repLatt.swap(mixer.repLatt);
        oldDensities.swap(mixer.oldDensities);
        residules.swap(mixer.residules);
        fft.swap(mixer.fft);
        std::swap(mixIteration, mixer.mixIteration);
    }

    template<class ScalarType, bool isSpinPolarized>
    void ChargeMixer<ScalarType, isSpinPolarized>::fetchInputDensity(DensityType& input) {
        oldDensities[mixIteration].swap(input);
    }

    template<class ScalarType, bool isSpinPolarized>
    void ChargeMixer<ScalarType, isSpinPolarized>::mix(const DensityType& output) {
        const DensityType& input = oldDensities[mixIteration];
        residules[mixIteration].getTotalDensity().flatten() = output.getTotalDensity().flatten() - input.getTotalDensity().flatten();
    }

    template<class ScalarType, bool isSpinPolarized>
    void ChargeMixer<ScalarType, isSpinPolarized>::mix(size_t iteration, DensityType& result) {
        if (iteration == 0) {
            using GridType = typename FFT3D::KSpaceType;
            using Index3D = typename GridType::Index3D;
            const auto& deltaRho = residules[0].getTotalDensity();
            fft.getRSpace().flatten() = deltaRho.flatten();
            fft.transform();

            auto& kSpace = fft.getKSpace();
            auto kernel = [&kSpace](Vector<ScalarType, 3> k, Index3D index) {
                const ScalarType kNorm = k.squaredNorm();
                const ScalarType factor = ScalarType(amix) * std::min(kNorm / (kNorm + square(ScalarType(bmix))), ScalarType(amin));
                kSpace(index) *= factor;
            };
            GridType::template forPointIndexInGrid<ScalarType, true, decltype(kernel)>(kSpace.getDim(), repLatt, kernel);
            fft.invTransform();

            const auto& rho_old = oldDensities[0].getTotalDensity().flatten();
            auto rho_new = result.getTotalDensity().flatten();
            rho_new = rho_old + fft.getRSpace().flatten();
        }
        else {
            using MatrixType = DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Element>;
            const size_t numValidRecord = iteration > mixIteration ? (DIISBufferSize - 1) : (mixIteration + 1);
            MatrixType diisMat = MatrixType(numValidRecord + 1, numValidRecord + 1, 1.0);
            diisMat(0, 0) = ScalarType(0);
            /* Construct equation */ {
                for (size_t i = 1; i < diisMat.getRow(); ++i) {
                    for (size_t j = i; j < diisMat.getCol(); ++j) {
                        ScalarType temp = residules[i - 1].getTotalDensity().flatten() * residules[j - 1].getTotalDensity().flatten();
                        diisMat(i, j) = temp;
                        diisMat(j, i) = temp;
                    }
                }
            }
            Vector<ScalarType> x(diisMat.getRow());
            /* Solve */ {
                auto b = Vector<ScalarType>(diisMat.getRow(), 0);
                b[0] = 1;
                const MatrixType inv_A = diisMat.inverse();
                x = inv_A * b;
            }

            auto rho_new = result.getTotalDensity().flatten();
            rho_new = ScalarType(0);
            for (size_t i = 1; i < x.getLength(); ++i) {
                const auto& rho_old = oldDensities[i - 1].getTotalDensity().flatten();
                const auto& residule = residules[i - 1].getTotalDensity().flatten();
                rho_new += (rho_old + ScalarType(pulay_mix) * residule) * x[i];
            }

            for (size_t i = 0; i < rho_new.getLength(); ++i)
                if (rho_new[i].isNegative())
                    rho_new[i] = 0;
        }
        mixIteration = (mixIteration + 1) % DIISBufferSize;
    }
}
