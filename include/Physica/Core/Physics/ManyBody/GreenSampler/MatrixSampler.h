/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "GreenSampler.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ImagKinetic.h"

namespace Physica {
    /**
     * Defination of AFM, refer to [1]
     *
     * Reference:
     * [1] Phys. Rev. Lett. 62, 591; https://doi.org/10.1103/PhysRevLett.62.591
     */
    template<Scalar T>
    class MatrixSampler : public GreenSampler<T> {
        using This = MatrixSampler<T>;
        using Base = GreenSampler<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;
        constexpr static int Dim = 2;

        using typename Base::Tr;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        enum Observable : char {
            Spin,
            Charge,
            DoubleOccupy,
        };
    private:
        LatticeModel<Dim> lattice;
        DenseTensor<T, 3> observes;
        Observable type;
    public:
        MatrixSampler(const HubbardParams<T>& params, const LatticeModel<Dim>& lattice, size_t numSample, Observable type);
        MatrixSampler(const This&) = delete;
        MatrixSampler(This&&) noexcept = delete;
        ~MatrixSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const GreenPair& greens, T rsign);

        [[nodiscard]] MatrixND<T> calcRawMean() const;
        [[nodiscard]] MatrixND<T> calcMean() const;
        [[nodiscard]] MatrixND<Tr> calcStructFactor() const;
        /* Getters */
        using Base::getNumSample;
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
        [[nodiscard]] size_t getNumSiteX() const noexcept { return lattice.getNumCellX(); }
        [[nodiscard]] size_t getNumSiteY() const noexcept { return lattice.getNumCellY(); }
        [[nodiscard]] size_t getNumSite() const noexcept { return getNumSiteX() * getNumSiteY(); }
    private:
        [[nodiscard]] MatrixND<T> calcCorrelation(const MatrixND<T>& greenU, const MatrixND<T>& greenD) const noexcept;
        [[nodiscard]] T calcCorrelation(const MatrixND<T>& greenU, int siteA, const MatrixND<T>& greenD, int siteB) const noexcept;

        using Base::calcDensityCorr;
    };

    template<Scalar T>
    MatrixSampler<T>::MatrixSampler(const HubbardParams<T>& params, const LatticeModel<Dim>& lattice, size_t numSample, Observable type)
            : Base(params, numSample)
            , lattice(lattice)
            , observes({numSample, getNumSiteX(), getNumSiteY()})
            , type(type) {}

    template<Scalar T>
    void MatrixSampler<T>::sample(const GreenPair& greens, T rsign) {
        observes.slice(Base::getCursor(), var(), var()) = calcCorrelation(greens[0], greens[1]) * rsign;
        Base::sample(rsign);
    }

    template<Scalar T>
    auto MatrixSampler<T>::calcRawMean() const -> MatrixND<T> {
        MatrixND<T> result(getNumSiteX(), getNumSiteY());
        for (size_t i = 0; i < observes.dim(0); ++i)
            result.toNextMean(i, observes.slice(Base::getCursor(), var(), var()));
        return result;
    }

    template<Scalar T>
    auto MatrixSampler<T>::calcMean() const -> MatrixND<T> {
        return calcRawMean() * reciprocal(Base::calcRSign());
    }

    template<Scalar T>
    auto MatrixSampler<T>::calcStructFactor() const -> MatrixND<Tr> {
        FFT<T, Dim> fft(lattice.getSuperSize(), PlanFlag::Estimate);
        fft.getRSpace() = calcMean();
        fft.transform();

        MatrixND<Tr> result = fft.getKSpace().squaredNorms();
        if (type == Spin) {
            // Add a minimum value to highlight the paramagnetic phase
            result += Trv(std::numeric_limits<T>::min());
        }
        else if (type == Charge) {
            // Ignore numerical errors
            const int kX = getNumSiteX();
            const int kY = FFT<T, 1>::rSizeToKSize(getNumSiteY());
            for (int x = 0; x < kX; ++x)
                for (int y = 0; y < kY; ++y)
                    if (result[x, y] < square(Trv(std::numeric_limits<T>::epsilon())))
                        result[x, y] = 0;
        }
        return result;
    }

    template<Scalar T>
    auto MatrixSampler<T>::calcCorrelation(const MatrixND<T>& greenU, const MatrixND<T>& greenD) const noexcept -> MatrixND<T> {
        MatrixND<T> result(getNumSiteX(), getNumSiteY());
        for (size_t siteA = 0; siteA < getNumSite(); ++siteA) {
            const auto indexA = lattice.toIndexND(siteA);
            for (size_t x = 0; x < getNumSiteX(); ++x) {
                for (size_t y = 0; y < getNumSiteY(); ++y) {
                    const auto indexB = indexA.shift(0, x, getNumSiteX()).shift(1, y, getNumSiteY());
                    const size_t siteB = lattice.toIndex1D(indexB);
                    result[x, y].toNextMean(siteA, calcCorrelation(greenU, siteA, greenD, siteB));
                }
            }
        }
        return result;
    }

    template<Scalar T>
    T MatrixSampler<T>::calcCorrelation(const MatrixND<T>& greenU, int siteA, const MatrixND<T>& greenD, int siteB) const noexcept {
        switch (type) {
        case Spin:
            return calcDensityCorr(greenU, siteA, siteB) + calcDensityCorr(greenD, siteA, siteB)
                 - calcDensityCorr(greenU, siteA, greenD, siteB) - calcDensityCorr(greenU, siteB, greenD, siteA);
        case Charge:
            return calcDensityCorr(greenU, siteA, siteB) + calcDensityCorr(greenD, siteA, siteB)
                 + calcDensityCorr(greenU, siteA, greenD, siteB) + calcDensityCorr(greenU, siteB, greenD, siteA);
        case DoubleOccupy:
            // TODO: Avoid mean field approximation
            return calcDensityCorr(greenU, siteA, greenD, siteB) * calcDensityCorr(greenU, siteB, greenD, siteA);
        default:
            unreachable();
        }
    }
}
