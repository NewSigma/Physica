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
        const MatrixND<T> hopUp = UnitMatrix<T>(getNumSite()) - greens[0];
        const MatrixND<T> hopDown = UnitMatrix<T>(getNumSite()) - greens[1];

        using IndexType = SquareLattice<2>::IndexType;
        auto buffer = MatrixND<T>(corr.getLength(), getNumSite());
        lattice.forSiteInLattice([&](const IndexType indexM) {
            const size_t dimX = lattice.getNumCellX();
            const size_t idM = lattice.toIndex1D(indexM);
            const size_t idN = lattice.toIndex1D(indexM.addX(1, dimX));
            for (size_t i = 0; i < dimX; ++i) {
                const size_t idI = lattice.toIndex1D(indexM.addX(i, dimX));
                const size_t idJ = lattice.toIndex1D(indexM.addX(i + 1, dimX));
                buffer(i, idM) = (hopUp(idM, idI) * hopDown(idN, idJ) + hopUp(idM, idJ) * hopDown(idN, idI)) * Trv(0.5);
            }
        });
        corr.toNextMean(Base::getCursor(), buffer.sum_cols() * reciprocal(Trv(getNumSite())));
        Base::sample(sign);
    }
}
