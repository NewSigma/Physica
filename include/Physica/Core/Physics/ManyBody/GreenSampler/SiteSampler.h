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

namespace Physica {
    template<Scalar T>
    class SiteSampler : public GreenSampler<T> {
        using This = SiteSampler<T>;
        using Base = GreenSampler<T>;
        using typename Base::Tv;
        using typename Base::Trv;

        using MatrixND = DQMC<T>::MatrixND;
        using GreenArray = DQMC<T>::GreenArray;
    public:
        enum Observable : char {
            Density,
            DoubleOccupy,
            Kinetic,
            Potential,
            Internal
        };
    private:
        VectorND<T> observes;
    public:
        SiteSampler(const DQMC<T>& dqmc_, size_t numSample);
        SiteSampler(const This&) = delete;
        SiteSampler(This&&) noexcept = delete;
        ~SiteSampler() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const Array<MatrixND, 2>& greens, Observable type);
        void sample(const GreenArray& greens, Observable type);

        [[nodiscard]] T calcMean() const;
        /* Getters */
        [[nodiscard]] const auto& getObserves() const noexcept { return observes; }
    private:
        T calcObservable(const MatrixND& greenU, const MatrixND& greenD, Observable type) const noexcept;
    };

    template<Scalar T>
    SiteSampler<T>::SiteSampler(const DQMC<T>& dqmc_, size_t numSample)
            : Base(dqmc_, numSample)
            , observes(numSample) {}

    template<Scalar T>
    void SiteSampler<T>::sample(const Array<MatrixND, 2>& greens, Observable type) {
        observes[Base::getCursor()] = calcObservable(greens[0], greens[1], type);
        Base::sample();
    }

    template<Scalar T>
    void SiteSampler<T>::sample(const GreenArray& greens, Observable type) {
        T mean = 0;
        for (int i = 0; i < greens.getCol(); ++i)
            mean.toNextMean(i, calcObservable(greens(0, i), greens(1, i), type));
        observes[Base::getCursor()] = mean;
        Base::sample();
    }

    template<Scalar T>
    T SiteSampler<T>::calcMean() const {
        return Base::calcMean(observes);
    }

    template<Scalar T>
    T SiteSampler<T>::calcObservable(const MatrixND& greenU, const MatrixND& greenD, Observable type) const noexcept {
        switch (type) {
        case Density:
            return T(2) - greenU.diag().reals().mean() - greenD.diag().reals().mean();
        case DoubleOccupy:
            return (T(1) - greenU.diag().reals()) * (T(1) - greenD.diag().reals()) / T(Base::getNumSite());
        case Kinetic: {
            const auto& hoppingT = Base::getHoppingMatrix();
            return -hadamard(hoppingT, greenU + greenD).sum().real() / T(Base::getNumSite());
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
