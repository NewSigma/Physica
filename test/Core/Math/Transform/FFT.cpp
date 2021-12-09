/*
 * Copyright 2021 WeiBo He.
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
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Utils/TestHelper.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;

ScalarType func(ScalarType x) {
    return sin(ScalarType(2 * M_PI * 3) * x) + sin(ScalarType(2 * M_PI * 4) * x) * 2;
}

int main() {
    const size_t N = 100;
    const double t_max = 2;
    
    Vector<ScalarType> data(N);
    {
        const Vector<ScalarType> x = Vector<ScalarType>::linspace(ScalarType::Zero(), ScalarType(t_max), N);
        for (size_t i = 0; i < N; ++i)
            data[i] = func(x[i]);
    }
    FFT<ScalarType> fft(data, ScalarType(t_max / N));
    fft.transform();
    const Vector<ScalarType> intense = toNormVector(fft.getComponents());

    /* Parseval theorem */ {
        const ScalarType power = square(data).sum();
        const ScalarType power_fft = square(intense).sum() / ScalarType(t_max * t_max / N);
        if (power != power_fft)
            return 1;
    }

    const ScalarType freq1_power = fft.getFreqIntense(3).norm();
    const ScalarType freq2_power = fft.getFreqIntense(4).norm();
    if (!scalarNear(freq2_power / freq1_power, ScalarType(2), 1E-14))
        return 1;
    return 0;
}
