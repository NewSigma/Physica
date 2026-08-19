/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ImagKinetic.h"

namespace Physica {
    /**
     * \class PairingSampler: Samples pairing-wave correlation function for square hubbard model
     *
     * Reference:
     * [1] Phys. Rev. B 39, 839(R) (1989); https://doi.org/10.1103/PhysRevB.39.839
     */
    template<Scalar T>
    class PairingSampler : public GreenSampler<T> {
        using This = PairingSampler<T>;
        using Base = GreenSampler<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;
        using IndexType = SquareLattice<2>::IndexType;
        using typename Base::Trv;
    public:
        // Refer to Fig.1 of [1] for pairing types
        enum Observable : char {
            Sxxyy,
            Dxxyy,
        };
    private:
        SquareLattice<2> lattice;
        MatrixND<T> corr;
        Observable type;
    public:
        PairingSampler(const HubbardParams<T>& params, SquareLattice<2> lattice, size_t numSample, Observable type);
        PairingSampler(const This&) = delete;
        PairingSampler(This&&) noexcept = delete;
        ~PairingSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const GreenPair& greens, T rsign);

        [[nodiscard]] MatrixND<T> calcRawMean() const;
        [[nodiscard]] MatrixND<T> calcMean() const;
        /* Getters */
        [[nodiscard]] size_t getNumSite() const noexcept { return lattice.getNumSuperCellSite(); }
    private:
        [[nodiscard]] T calcWavePairingCorr(const GreenPair& greens, IndexType from, IndexType to);
        [[nodiscard]] VectorND<T> makeWavePairingWeight();
        void forWaveNeighbor(IndexType center, std::invocable<int, IndexType> auto fn);
        /* Static members */
        [[nodiscard]] static T calcPairingCorr(const GreenPair& greens, int siteM, int siteN, int siteI, int siteJ);
    };

    template<Scalar T>
    PairingSampler<T>::PairingSampler(const HubbardParams<T>& params, SquareLattice<2> lattice, size_t numSample, Observable type)
            : Base(params, numSample), lattice(std::move(lattice)), corr(lattice.getNumCellX(), lattice.getNumCellY()), type(type) {}

    template<Scalar T>
    void PairingSampler<T>::sample(const GreenPair& greens, T rsign) {
        auto buffer = MatrixND<T>(corr.getRow(), corr.getCol(), 0);
        for (int siteI = 0; siteI < getNumSite(); ++siteI) {
            const auto from = lattice.toIndexND(siteI);
            for (int siteM = 0; siteM < getNumSite(); ++siteM) {
                const auto to = lattice.toIndexND(siteM);
                auto dist = to.shift(0, -ssize_t(from[0]), lattice.getNumCellX())
                              .shift(1, -ssize_t(from[1]), lattice.getNumCellY());
                buffer[dist.getX(), dist.getY()] += calcWavePairingCorr(greens, from, to) / Trv(getNumSite());
            }
        }
        corr.toNextMean(Base::getCursor(), buffer);
        Base::sample(rsign);
    }

    template<Scalar T>
    MatrixND<T> PairingSampler<T>::calcRawMean() const {
        return corr;
    }

    template<Scalar T>
    MatrixND<T> PairingSampler<T>::calcMean() const {
        return corr * reciprocal(Base::calcRSign());
    }

    template<Scalar T>
    T PairingSampler<T>::calcWavePairingCorr(const GreenPair& greens, IndexType from, IndexType to) {
        auto weight = makeWavePairingWeight();
        MatrixND<T> kernel(weight.getLength());
        forWaveNeighbor(from, [&](int r, IndexType fromNeigh) {
            forWaveNeighbor(to, [&](int c, IndexType toNeigh) {
                const int siteI = lattice.toIndex1D(from);
                const int siteJ = lattice.toIndex1D(fromNeigh);
                const int siteM = lattice.toIndex1D(to);
                const int siteN = lattice.toIndex1D(toNeigh);
                kernel[r, c] = calcPairingCorr(greens, siteM, siteN, siteI, siteJ);
            });
        });
        return weight.transpose() * kernel * weight / T(kernel.getSize());
    }

    template<Scalar T>
    VectorND<T> PairingSampler<T>::makeWavePairingWeight() {
        switch (type) {
        case Sxxyy:
            return {1, 1, 1, 1};
        case Dxxyy:
            return {1, -1, 1, -1};
        default:
            unreachable();
        }
    }

    template<Scalar T>
    void PairingSampler<T>::forWaveNeighbor(IndexType center, std::invocable<int, IndexType> auto fn) {
        switch (type) {
        case Sxxyy:
        case Dxxyy: {
            Array<int, 4> shiftDims{0, 1, 0, 1};
            Array<int, 4> deltas{1, 1, -1, -1};
            for (int i = 0; i < shiftDims.size(); ++i) {
                int shiftDim = shiftDims[i];
                fn(i, center.shift(shiftDim, deltas[i], lattice.getSuperSize()[shiftDim]));
            }
            break;
        }
        default:
            unreachable();
        }
    }

    template<Scalar T>
    T PairingSampler<T>::calcPairingCorr(const GreenPair& greens, int siteM, int siteN, int siteI, int siteJ) {
        auto calcHopping = [](const MatrixND<T>& green, int siteA, int siteB) static noexcept {
            return Trv(siteA == siteB) - green[siteB, siteA];
        };

        const auto& [greenU, greenD] = greens;
        T x1 = fma(calcHopping(greenU, siteM, siteI), calcHopping(greenD, siteN, siteJ), Trv(0));
        T x2 = fma(calcHopping(greenU, siteM, siteJ), calcHopping(greenD, siteN, siteI), x1);
        T x3 = fma(calcHopping(greenD, siteM, siteI), calcHopping(greenU, siteN, siteJ), x2);
        return fma(calcHopping(greenD, siteM, siteJ), calcHopping(greenU, siteN, siteI), x3);
    }
}
