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
