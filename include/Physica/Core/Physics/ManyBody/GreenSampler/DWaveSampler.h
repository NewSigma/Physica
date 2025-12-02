/*
 * Copyright 2025 Weibo He.
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

#include "GreenSampler.h"
#include "Physica/Core/Math/Statistics/Correlation.h"
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ImagKinetic.h"

namespace Physica {
    /**
     * \class DWaveSampler: Samples d-wave pairing correlation function for square hubbard model
     */
    template<Scalar T>
    class DWaveSampler : public GreenSampler<T> {
        using This = DWaveSampler<T>;
        using Base = GreenSampler<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;
        using typename Base::Trv;

        SquareLattice<2> lattice;
        VectorND<T> corr;
    public:
        DWaveSampler(const SquareLattice<2>& lattice, const HubbardParams<T>& params, size_t numSample);
        DWaveSampler(const This&) = delete;
        DWaveSampler(This&&) noexcept = delete;
        ~DWaveSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const GreenPair& greens, T sign);
        /* Getters */
        [[nodiscard]] size_t getNumSite() const noexcept { return lattice.getNumSuperCellSite(); }
        [[nodiscard]] const auto& getCorr() const noexcept { return corr; }
    };

    template<Scalar T>
    DWaveSampler<T>::DWaveSampler(const SquareLattice<2>& lattice, const HubbardParams<T>& params, size_t numSample)
            : Base(params, numSample), lattice(lattice), corr(lattice.getNumCellX()) {}

    template<Scalar T>
    void DWaveSampler<T>::sample(const GreenPair& greens, T sign) {
        auto buffer = MatrixND<T>(corr.getLength(), getNumSite());
        for (int siteM = 0; siteM < getNumSite(); ++siteM) {
            const auto indexM = lattice.toIndexND(siteM);
            const size_t dimX = lattice.getNumCellX();
            const size_t siteN = lattice.toIndex1D(indexM.shift(0, 1, dimX));
            for (size_t i = 0; i < dimX; ++i) {
                using Base::calcDensityCorr;
                const size_t siteI = lattice.toIndex1D(indexM.shift(0, i, dimX));
                const size_t siteJ = lattice.toIndex1D(indexM.shift(0, i + 1, dimX));
                buffer(i, siteM) = calcDensityCorr(greens[0], siteM, siteI) * calcDensityCorr(greens[1], siteN, siteJ)
                                 + calcDensityCorr(greens[0], siteM, siteJ) * calcDensityCorr(greens[1], siteN, siteI)
                                 + calcDensityCorr(greens[1], siteM, siteI) * calcDensityCorr(greens[0], siteN, siteJ)
                                 + calcDensityCorr(greens[1], siteM, siteJ) * calcDensityCorr(greens[0], siteN, siteI);
            }
        }
        buffer *= Trv(0.5);
        corr.toNextMean(Base::getCursor(), buffer.sum_cols() * reciprocal(Trv(getNumSite())));
        Base::sample(sign);
    }
}
