/*
 * Copyright 2021-2022 WeiBo He.
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
#include <iostream>
#include "Physica/Core/Math/Transform/FFT.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using RealType = Scalar<Double, false>;
using ComplexType = ComplexScalar<RealType>;

int main() {
    /* 1D real */ {
        const size_t N = 100;
        const double t_max = 2;
        constexpr double freq1 = 3;
        constexpr double freq2 = 4;
        
        Vector<RealType> data(N);
        {
            const Vector<RealType> v_x = Vector<RealType>::linspace(RealType::Zero(), RealType(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(RealType(2 * M_PI * freq1) * x) + sin(RealType(2 * M_PI * freq2) * x) * 2;
            }
        }
        FFT<RealType> fft(data, RealType(t_max / N), FFT<RealType>::Measure);
        const Vector<RealType> intense = toNormVector(fft.getKSpace());

        /* Parseval theorem */ {
            const RealType power = square(data).sum();
            const RealType power_fft = square(intense).sum() / RealType(intense.getLength() - 1);
            if (!scalarNear(power, power_fft, 1E-15))
                return 1;
        }
        /* Test freq */ {
            const double kSpaceDelta = double(fft.getKSpaceDelta());
            const RealType freq1_power = intense[freq1 / kSpaceDelta];
            const RealType freq2_power = intense[freq2 / kSpaceDelta];
            if (!scalarNear(freq2_power / freq1_power, RealType(2), 1E-14))
                return 1;
        }
        /* Test inv */ {
            constexpr double precision = 1E-13;
            fft.invTransform(fft.getKSpace());
            for (size_t i = 0; i < data.getLength(); ++i) {
                const bool isNear = scalarNear(data[i], fft.getRSpace()[i], 1E-14);
                const bool isSmall = abs(data[i]) < RealType(precision) && abs(fft.getRSpace()[i]) < RealType(precision);
                if(!isNear && !isSmall)
                    return 1;
            }
        }
    }
    /* 1d complex */ {
        const size_t N = 100;
        const double t_max = 2;
        constexpr double freq1 = 3;
        constexpr double freq2 = 4;
        
        Vector<ComplexType> data(N);
        {
            const Vector<RealType> v_x = Vector<RealType>::linspace(RealType::Zero(), RealType(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(RealType(2 * M_PI * freq1) * x) + sin(RealType(2 * M_PI * freq2) * x) * 2;
            }
        }
        FFT<ComplexType> fft(data, RealType(t_max / N), FFT<ComplexType>::Estimate);
        Vector<ComplexType> trans(N);
        for (size_t i = 0; i < N; ++i) {
            ComplexType temp(0);
            for (size_t j = 0; j < N; ++j) {
                RealType phase = 2 * M_PI * j * i / N;
                temp += data[j] * ComplexType(cos(phase), -sin(phase));
            }
            trans[i] = temp;
        }
        fft.invTransform(fft.getKSpace());
        if (!vectorNear(data, fft.getRSpace(), 1E-14))
            return 1;
    }
    /* 2d real */ {
        const size_t N1 = 50;
        const size_t N2 = 100;
        const double deltaX = 0.01;
        const double deltaY = 0.01;
        constexpr double freq1 = 10;
        constexpr double freq2 = 5;

        Vector<RealType> data(N1 * N2);
        {
            size_t index = 0;
            for (size_t i = 0; i < N1; ++i) {
                for (size_t j = 0; j < N2; ++j) {
                    data[index++] = RealType(std::sin(2 * M_PI * freq1 * i * deltaX) + 2 * std::cos(2 * M_PI * freq2 * j * deltaY));
                }
            }
        }
        FFT<RealType, 2> fft(data, {N1, N2}, {deltaX, deltaY});
        /* Test freq */ {
            auto intense = toNormVector(fft.getKSpace());
            const RealType freq1_power = intense[size_t(freq1 / double(fft.getKSpaceDelta(0))) * fft.getKSpaceSize(1)];
            const RealType freq1_power_conj = intense[(N1 - size_t(freq1 / double(fft.getKSpaceDelta(0)))) * fft.getKSpaceSize(1)];
            const RealType freq2_power = intense[freq2 / double(fft.getKSpaceDelta(1))];
            if (!scalarNear(freq1_power, freq1_power_conj, 1E-15))
                return 1;
            if (!scalarNear(freq2_power / freq1_power, RealType(2), 1E-14))
                return 1;
        }
        /* Test inv */ {
            constexpr double precision = 1E-11;
            fft.invTransform(fft.getKSpace());
            for (size_t i = 0; i < data.getLength(); ++i) {
                const bool isNear = scalarNear(data[i], fft.getRSpace()[i], precision);
                const bool isSmall = abs(data[i]) < RealType(precision) && abs(fft.getRSpace()[i]) < RealType(precision);
                if(!isNear && !isSmall)
                    return 1;
            }
        }
    }
    return 0;
}
