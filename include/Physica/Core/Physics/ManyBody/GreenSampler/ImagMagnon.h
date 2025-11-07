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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"

namespace Physica {
    template<Scalar T>
    class ImagMagnon : public GreenSampler<T> {
        using This = ImagMagnon<T>;
        using Base = GreenSampler<T>;
        using Tc = T::ComplexType;
    private:
        DenseTensor<Tc, 3> magnons; 
        FFT<T> fft;
    public:
        ImagMagnon(const DQMC<T>& dqmc, size_t numSample);
        ImagMagnon(const This&) = delete;
        ImagMagnon(This&&) noexcept = delete;
        ~ImagMagnon() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void sample(const MatrixND<T>& aux);
        [[nodiscard]] MatrixND<T> calcMean() const;
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return magnons.dim(1); }
        [[nodiscard]] int getNumSplit() const noexcept { return fft.getRSpaceSize(); }
    };

    template<Scalar T>
    ImagMagnon<T>::ImagMagnon(const DQMC<T>& dqmc, size_t numSample)
            : Base(dqmc, numSample)
            , magnons(numSample, dqmc.getNumSite(), FFT<T>::rSizeToKSize(dqmc.getNumSplit()))
            , fft(dqmc.getNumSplit(), PlanFlag::Estimate) {}

    template<Scalar T>
    void ImagMagnon<T>::sample(const MatrixND<T>& aux) {
        for (size_t site = 0; site < aux.getRow(); ++site) {
            fft.transform(aux.row(site));
            magnons.fiber(2, {Base::getCursor(), site, Dynamic}) = fft.getKSpace() * reciprocal(T(getNumSplit()));
        }
        Base::sample();
    }

    template<Scalar T>
    auto ImagMagnon<T>::calcMean() const -> MatrixND<T> {
        MatrixND<T> result(getNumSite(), magnons.dim(2));
        for (size_t r = 0; r < result.getRow(); ++r)
            for (size_t c = 0; c < result.getCol(); ++c)
                result(r, c) = magnons.fiber(0, {Dynamic, r, c}).squaredNorms().mean();
        result *= reciprocal(Base::calcSign());
        return result;
    }
}
