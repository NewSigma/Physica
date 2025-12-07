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
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ImagKinetic.h"

namespace Physica {
    template<Scalar T>
    class SiteSampler : public GreenSampler<T> {
        using This = SiteSampler<T>;
        using Base = GreenSampler<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;

        using typename Base::Tv;
        using typename Base::Trv;
    public:
        enum Observable : char {
            Density,
            Density2,
            DoubleOccupy,
            MagMoment,
            Kinetic,
            Potential,
            Internal,
        };
    private:
        VectorND<T> observes;
        Observable type;
    public:
        SiteSampler(const HubbardParams<T>& params, size_t numSample, Observable type);
        SiteSampler(const This&) = delete;
        SiteSampler(This&&) noexcept = delete;
        ~SiteSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const GreenPair& greens, T rsign);

        [[nodiscard]] T calcMean() const;
        /* Getters */
        using Base::getNumSite;
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
    private:
        [[nodiscard]] T calcObservable(const MatrixND<T>& greenU, const MatrixND<T>& greenD, Observable typeX) const noexcept;
    };

    template<Scalar T>
    SiteSampler<T>::SiteSampler(const HubbardParams<T>& params, size_t numSample, Observable type)
            : Base(params, numSample)
            , observes(numSample)
            , type(type) {}

    template<Scalar T>
    void SiteSampler<T>::sample(const GreenPair& greens, T rsign) {
        observes[Base::getCursor()] = calcObservable(greens[0], greens[1], type) * rsign;
        Base::sample(rsign);
    }

    template<Scalar T>
    T SiteSampler<T>::calcMean() const {
        return observes.mean() / Base::calcSign();
    }

    template<Scalar T>
    T SiteSampler<T>::calcObservable(const MatrixND<T>& greenU, const MatrixND<T>& greenD, Observable typeX) const noexcept {
        switch (typeX) {
        case Density:
            return T(2) - greenU.diag().reals().mean() - greenD.diag().reals().mean();
        case Density2: {
            int numSite = getNumSite();
            T result = 0;
            for (int siteA = 0; siteA < numSite; ++siteA) {
                for (int siteB = 0; siteB < numSite; ++siteB) {
                    result += Base::calcDensityCorr(greenU, siteA, siteB);
                    result += Base::calcDensityCorr(greenD, siteA, siteB);
                }
            }
            result /= Trv(numSite * numSite);

            T rhoU = T(1) - greenU.diag().reals().mean();
            T rhoD = T(1) - greenD.diag().reals().mean();
            result += rhoU * rhoD * Trv(2);
            return result;
        }
        case DoubleOccupy:
            return (T(1) - greenU.diag().reals()) * (T(1) - greenD.diag().reals()) / T(getNumSite());
        case MagMoment:
            return calcObservable(greenU, greenD, Density) - calcObservable(greenU, greenD, DoubleOccupy) * Trv(2);
        case Kinetic: {
            const auto& hoppingT = Base::getHoppingMatrix();
            return -hadamard(hoppingT, greenU + greenD).sum().real() / T(getNumSite());
        }
        case Potential:
            return calcObservable(greenU, greenD, DoubleOccupy) * Base::getRepelU();
        case Internal:
            return calcObservable(greenU, greenD, Kinetic) + calcObservable(greenU, greenD, Potential);
        default:
            unreachable();
        }
    }
}
