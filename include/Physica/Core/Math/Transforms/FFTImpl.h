/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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
    FFT<type, errorTrack>::FFT(const FFT& fft) : data(fft.data) {}

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>::FFT(FFT&& fft) noexcept : data(std::move(fft.data)) {}

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>& FFT<type, errorTrack>::operator=(const FFT& fft) {
        data = fft.data;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    FFT<type, errorTrack>& FFT<type, errorTrack>::operator=(FFT&& fft) noexcept {
        data = std::move(fft.data);
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    inline void FFT<type, errorTrack>::transform() {
        transformImpl((MathConst::getInstance() / data.getLength()) << 1);
    }

    template<ScalarType type, bool errorTrack>
    inline void FFT<type, errorTrack>::invTransform() {
        transformImpl(-(MathConst::getInstance() / data.getLength()) << 1);
    }

    template<ScalarType type, bool errorTrack>
    void FFT<type, errorTrack>::transformImpl(const Scalar<type, errorTrack>&& phase) {
        data *= distance;
        const auto length = data.getLength();
        CStyleArray<ComplexScalar<type, errorTrack>, Dynamic> array(length);
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
        data = std::move(array);
    }
}

#endif