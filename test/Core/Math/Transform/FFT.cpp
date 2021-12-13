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
using RealType = Scalar<Double, false>;
using ComplexType = ComplexScalar<RealType>;

int main() {
    /* 1D real */ {
        const size_t N = 100;
        const double t_max = 2;
        
        Vector<RealType> data(N);
        {
            const Vector<RealType> v_x = Vector<RealType>::linspace(RealType::Zero(), RealType(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(RealType(2 * M_PI * 3) * x) + sin(RealType(2 * M_PI * 4) * x) * 2;
            }
        }
        FFT<RealType> fft(data, RealType(t_max / N));
        fft.transform();
        const Vector<RealType> intense = toNormVector(fft.getComponents());

        /* Parseval theorem */ {
            const RealType power = square(data).sum();
            const RealType power_fft = square(intense).sum() / RealType(fft.getDeltaT() * fft.getDeltaT() * (intense.getLength() - 1));
            if (!scalarNear(power, power_fft, 1E-15))
                return 1;
        }
        const RealType freq1_power = fft.getFreqIntense(3).norm();
        const RealType freq2_power = fft.getFreqIntense(4).norm();
        if (!scalarNear(freq2_power / freq1_power, RealType(2), 1E-14))
            return 1;
    }
    /* 2d real */ {
        const size_t N1 = 50;
        const size_t N2 = 100;
        const double deltaX = 0.01;
        const double deltaY = 0.01;

        Vector<RealType> data(N1 * N2);
        {
            size_t index = 0;
            for (size_t i = 0; i < N1; ++i) {
                for (size_t j = 0; j < N2; ++j) {
                    data[index++] = RealType(std::sin(2 * M_PI * 10 * i * deltaX) + 2 * std::cos(2 * M_PI * 5 * j * deltaY));
                }
            }
        }
        FFT<RealType, 2> fft(data, {N1, N2}, {deltaX, deltaY});
        fft.transform();

        const RealType freq1_power = fft.getFreqIntense({10, 0}).norm();
        const RealType freq2_power = fft.getFreqIntense({0, 5}).norm();
        if (!scalarNear(freq2_power / freq1_power, RealType(2), 2E-2))
            return 1;
    }
    return 0;
}
