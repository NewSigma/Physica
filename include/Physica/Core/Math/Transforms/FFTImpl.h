/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_FFTIMPL_H
#define PHYSICA_FFTIMPL_H
/*!
 * This file is part of implementations of \FFT.
 * Do not include this header file, include FFT.h instead.
 */
namespace Physica::Core {
    /*!
     * Improve:
     * 1.Integrate methods may be various.
     * 2.Sampling distance does not have to be equal.
     * 3.A lot of other features.
     */
    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>::FFT(Vector<ComplexScalar<type, errorTrack>, Dynamic> data
            , Scalar<type, false> distance) : data(std::move(data)), distance(std::move(distance)) {}

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>::FFT(const FFT& fft) : data(fft.data) {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
    }

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>::FFT(FFT&& fft) noexcept : data(std::move(fft.data)) {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
    }

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>& FFT<type, errorTrack>::operator=(const FFT& fft) {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        data = fft.data;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>& FFT<type, errorTrack>::operator=(FFT&& fft) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        data = std::move(fft.data);
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    inline void FFT<type, errorTrack>::transform() {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        transformImpl((MathConst::getInstance() / data.getLength()) << 1);
    }

    template<ScalarType type, bool errorTrack>
    inline void FFT<type, errorTrack>::invTransform() {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        transformImpl(-(MathConst::getInstance() / data.getLength()) << 1);
    }

    template<ScalarType type, bool errorTrack>
    void FFT<type, errorTrack>::transformImpl(const Scalar<type, errorTrack>&& phase) {
        data *= distance;
        const auto length = data.getLength();
        CStyleArray<ComplexScalar<type, errorTrack>> array(length);
        //Optimize:
        //1.i and j is changeable.(dynamic programming)
        //2.Use the formula such as sin(a + b) to avoid calculate sin and cos directly.
        for(size_t i = 0; i < length; ++i) {
            const auto phase1 = phase * i;
            auto result_i = ComplexScalar<type, errorTrack>::getZero();
            for(size_t j = 0; j < length; ++j) {
                const auto phase2 = phase1 * j;
                result_i += ComplexScalar<type, errorTrack>(cos(phase2), sin(phase2)) * data[j];
            }
            array.allocate(std::move(result_i), i);
        }
        array.setLength(length);
        data = std::move(array);
    }
}

#endif