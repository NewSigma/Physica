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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"

namespace Physica {
    template<Scalar T>
    class ImagMagnon : public GreenSampler<T> {
        using This = ImagMagnon<T>;
        using Base = GreenSampler<T>;
        using Tc = T::ComplexType;
        using Trv = Base::Trv;
    private:
        DenseTensor<Tc, 3> magnons;
        FFT<T> fft;
    public:
        ImagMagnon(const HubbardParams<T>& params, size_t numSample);
        ImagMagnon(const This&) = delete;
        ImagMagnon(This&&) noexcept = delete;
        ~ImagMagnon() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const MatrixND<T>& aux, Trv rsign);
        [[nodiscard]] MatrixND<T> calcMean() const;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return magnons.dim(1); }
        [[nodiscard]] int getNumSplit() const noexcept { return fft.getRSpaceSize(); }
    };

    template<Scalar T>
    ImagMagnon<T>::ImagMagnon(const HubbardParams<T>& params, size_t numSample)
            : Base(params, numSample)
            , magnons(numSample, params.getNumSite(), FFT<T>::rSizeToKSize(params.getNumSplit()))
            , fft(params.getNumSplit(), PlanFlag::Estimate) {}

    template<Scalar T>
    void ImagMagnon<T>::sample(const MatrixND<T>& aux, Trv rsign) {
        for (size_t site = 0; site < aux.getRow(); ++site) {
            fft.transform(aux.row(site));
            magnons.fiber(Base::getCursor(), site, var()) = fft.getKSpace() * reciprocal(sqrt(T(getNumSplit())));
        }
        Base::sample(rsign);
    }

    template<Scalar T>
    auto ImagMagnon<T>::calcMean() const -> MatrixND<T> {
        auto result = MatrixND<T>::generate([&](size_t r, size_t c) {
            return magnons.fiber(var(), r, c).squaredNorms().mean();
        }, getNumSite(), magnons.dim(2));
        result *= reciprocal(Base::calcRSign());
        return result;
    }
}
