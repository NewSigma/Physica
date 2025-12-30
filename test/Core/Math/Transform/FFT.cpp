/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Transform/DiffFFT.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = T::ComplexType;

namespace {
    void testReal2D() {
        const size_t N1 = 50;
        const size_t N2 = 100;
        const double deltaX = 0.01;
        const double deltaY = 0.01;
        constexpr T freq1 = 10;
        constexpr T freq2 = 5;

        DenseMatrix<T, MatrixOption::Row> data(N1, N2);
        for (size_t i = 0; i < N1; ++i)
            for (size_t j = 0; j < N2; ++j)
                data[i, j] = T(sin(T(2 * i) * MathConst<T>::pi * freq1 * deltaX) + cos(T(2 * j) * MathConst<T>::pi * freq2 * deltaY) * 2);

        FFT<T, 2> fft({N1, N2}, PlanFlag::Estimate);
        fft.transform(data);
        /* Test freq */ {
            const T deltaFreq1 = fft.getKSpaceDelta(deltaX, 0) / (MathConst<T>::pi * 2);
            const T deltaFreq2 = fft.getKSpaceDelta(deltaY, 1) / (MathConst<T>::pi * 2);
            const auto index1 = size_t((freq1 / deltaFreq1).toMachine());
            const auto index2 = size_t((freq2 / deltaFreq2).toMachine());
            const VectorND<T> intense = fft.getKSpace().flatten().norms();
            const T freq1_power = intense[index1 * fft.getKSpaceSize()[1]];
            const T freq1_power_conj = intense[(N1 - index1) * fft.getKSpaceSize()[1]];
            const T freq2_power = intense[index2];
            expect(scalarNear(freq1_power, freq1_power_conj, 1E-15));
            expect(scalarNear(freq2_power / freq1_power, T(2), 1E-14));
        }
        /* Test inv */ {
            constexpr double precision = 1E-11;
            fft.invTransform();
            for (size_t i = 0; i < data.getRow(); ++i) {
                for (size_t j = 0; j < data.getCol(); ++j) {
                    const bool isNear = scalarNear(data[i, j], fft.getRSpace()[i, j], precision);
                    const bool isSmall = abs(data[i, j]) < T(precision) && abs(fft.getRSpace()[i, j]) < T(precision);
                    expect(isNear || isSmall);
                }
            }
        }
    }

    void test_differentiable() {
        using dfloat = Diff<T, DiffMode::Forward, 1>;
        using dcfloat = Diff<Tc, DiffMode::Forward, 1>;
        const size_t N = 100;
        constexpr double freq1 = 3;
        constexpr double freq2 = 4;

        VectorND<T> values(N);
        VectorND<T> grads(N);
        VectorND<dfloat> data(N);
        for (size_t i = 0; i < N; ++i) {
            const T x = T(i) * 0.01;
            values[i] = sin(T(2 * M_PI * freq1) * x) + sin(T(2 * M_PI * freq2) * x) * 2;
            grads[i] = cos(T(2 * M_PI * freq1) * x) * 2 + cos(T(2 * M_PI * freq2) * x);
            data[i] = dfloat(values[i], grads[i]);
        }
        VectorND<dcfloat> answer{};
        /* Make answer */ {
            FFT<T> fft(N, PlanFlag::Estimate);
            fft.getRSpace() = values;
            fft.transform();
            VectorND<Tc> k_values = fft.getKSpace();
            fft.getRSpace() = grads;
            fft.transform();
            VectorND<Tc> k_grads = fft.getKSpace();
            expect(k_values.getLength() == k_grads.getLength());

            answer.resize(k_values.getLength());
            for (size_t i = 0; i < answer.getLength(); ++i)
                answer[i] = dcfloat(k_values[i], k_grads[i]);
        }
        FFT<dfloat> fft(N, PlanFlag::Estimate);
        /* Test transform */ {
            fft.getRSpace() = data;
            fft.transform();
            expect(vectorNear(fft.getKSpace(), answer, 1E-15));
        }
        /* Test invTransform */ {
            constexpr double precision = 1E-13;
            fft.invTransform();
            for (size_t i = 0; i < data.getLength(); ++i) {
                const bool isNear = scalarNear(dfloat(data[i]), dfloat(fft.getRSpace()[i]), precision);
                const bool isSmall = abs(data[i].value()) < T(precision) && abs(fft.getRSpace()[i].value()) < T(precision);
                expect(isNear || isSmall);
            }
        }
    }
}

int main() {
    /* 1D real */ {
        const size_t N = 100;
        const double t_max = 2;
        constexpr double freq1 = 3;
        constexpr double freq2 = 4;
        
        VectorND<T> data(N);
        {
            const VectorND<T> v_x = VectorND<T>::linspace(T(0), T(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(T(2 * M_PI * freq1) * x) + sin(T(2 * M_PI * freq2) * x) * 2;
            }
        }
        FFT<T> fft(data, PlanFlag::Measure);
        const VectorND<T> intense = fft.getKSpace().norms();

        /* Parseval theorem */ {
            const T energyR = data.squaredNorm();
            const T energyK = fft.getKSpace().parseval();
            expect(scalarNear(energyR, energyK, 1E-15));
        }
        /* Test freq */ {
            const double deltaFreq = double(fft.getKSpaceDelta(T(t_max / N))) / (2 * M_PI);
            const T freq1_power = intense[size_t(freq1 / deltaFreq)];
            const T freq2_power = intense[size_t(freq2 / deltaFreq)];
            expect(scalarNear(freq2_power / freq1_power, T(2), 1E-14));
        }
        /* Test inv */ {
            constexpr double precision = 1E-13;
            fft.invTransform();
            for (size_t i = 0; i < data.getLength(); ++i) {
                const bool isNear = scalarNear(data[i], fft.getRSpace()[i], 1E-14);
                const bool isSmall = abs(data[i]) < T(precision) && abs(fft.getRSpace()[i]) < T(precision);
                expect(isNear || isSmall);
            }
        }
    }
    /* 1d complex */ {
        const size_t N = 100;
        const double t_max = 2;
        constexpr double freq1 = 3;
        constexpr double freq2 = 4;
        
        VectorND<Tc> data(N);
        {
            const VectorND<T> v_x = VectorND<T>::linspace(T(0), T(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(T(2 * M_PI * freq1) * x) + sin(T(2 * M_PI * freq2) * x) * 2;
            }
        }
        FFT<Tc> fft(data, PlanFlag::Estimate);
        VectorND<Tc> trans(N);
        for (size_t i = 0; i < N; ++i) {
            Tc temp(0);
            for (size_t j = 0; j < N; ++j) {
                T phase = MathConst<T>::pi * 2 * j * i / N;
                temp += data[j] * Tc(cos(phase), -sin(phase));
            }
            trans[i] = temp;
        }
        fft.invTransform();
        expect(vectorNear(data, fft.getRSpace(), 1E-14));
    }
    testReal2D();
    test_differentiable();
    return 0;
}
