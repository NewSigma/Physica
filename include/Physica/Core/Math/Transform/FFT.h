/*
 * Copyright 2020-2021 WeiBo He.
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

#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType> class FFTImpl;
    }

    template<class ScalarType>
    class FFT {
        using Impl = Internal::FFTImpl<ScalarType>;
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
    private:
        Impl impl;
    public:
        FFT(const Vector<ScalarType>& data, const RealType& deltaT);
        FFT(const FFT&) = default;
        FFT(FFT&&) noexcept = default;
        ~FFT() = default;
        /* Operators */
        FFT& operator=(FFT fft);
        ComplexType operator()(size_t i) { return impl(i); }
        /* Operations */
        void transform() { impl.transform(); }
        void invTransform() { impl.invTransform(); }
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return impl.getSize(); }
        [[nodiscard]] const RealType& getDeltaT() const noexcept { return impl.getDeltaT(); }
        [[nodiscard]] ComplexType getComponent(ssize_t index) const { return impl.getComponent(index); }
        [[nodiscard]] Vector<ComplexType> getComponents() const { return impl.getComponents(); }
        [[nodiscard]] RealType getDeltaFreq() const noexcept { return reciprocal(impl.getDistance() * impl.getSize()); }
        [[nodiscard]] ComplexType getFreqIntense(const RealType& freq) const noexcept;
        /* Helpers */
        void swap(FFT& fft) { impl.swap(fft.impl); }
    };

    template<class ScalarType>
    FFT<ScalarType>::FFT(const Vector<ScalarType>& data, const RealType& deltaT) : impl(data, deltaT) {}

    template<class ScalarType>
    FFT<ScalarType>& FFT<ScalarType>::operator=(FFT<ScalarType> fft) {
        swap(fft);
        return *this;
    }

    template<class ScalarType>
    typename FFT<ScalarType>::ComplexType FFT<ScalarType>::getFreqIntense(const RealType& freq) const noexcept {
        const double round_helper = freq.isPositive() ? 0.5 : -0.5;
        const double float_index = double(getDeltaT() * freq * RealType(getSize()));
        const ssize_t index = static_cast<ssize_t>(float_index + round_helper);
        return getComponent(index);
    }
}

#ifdef PHYSICA_FFTW3
    #include "FFTImpl_FFTW3.h"
#else
    #include "FFTImpl.h"
#endif
