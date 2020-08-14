/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FFT_H
#define PHYSICA_FFT_H

#include "Physica/Core/MultiPrecition/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"

namespace Physica::Core {
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class FFT {
    protected:
        Vector<ComplexScalar<type, errorTrack>, Dynamic> data;
        Scalar<type, false> distance;
    public:
        FFT(Vector<ComplexScalar<type, errorTrack>, Dynamic> data, Scalar<type, false> distance);
        FFT(const FFT& fft);
        FFT(FFT&& fft) noexcept;
        ~FFT() = default;
        /* Operators */
        FFT& operator=(const FFT& fft);
        FFT& operator=(FFT&& fft) noexcept;
        ComplexScalar<type, errorTrack> operator()(size_t i) { return data[i]; }
        /* Transforms */
        inline void transform();
        inline void invTransform();
    private:
        void transformImpl(const Scalar<type, errorTrack>&& phase);
    };
}

#include "FFTImpl.h"

#endif
