/*
 * Copyright 2021-2024 Weibo He.
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
#include "Physica/Core/Math/Transform/DiffFFT.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;
using RealType = float64;
using ComplexType = Complex<RealType>;

void test_differentiable() {
    using ValueType = float64;
    using ComplexPlainScalar = Complex<ValueType>;
    using ScalarType = Diff<ValueType, DiffMode::Forward, 1>;
    using ComplexType = Diff<ComplexPlainScalar, DiffMode::Forward, 1>;
    const size_t N = 100;
    constexpr double freq1 = 3;
    constexpr double freq2 = 4;

    VectorND<ValueType> values(N);
    VectorND<ValueType> grads(N);
    VectorND<ScalarType> data(N);
    for (size_t i = 0; i < N; ++i) {
        const ValueType x = ValueType(i) * 0.01;
        values[i] = sin(ValueType(2 * M_PI * freq1) * x) + sin(ValueType(2 * M_PI * freq2) * x) * 2;
        grads[i] = cos(ValueType(2 * M_PI * freq1) * x) * 2 + cos(ValueType(2 * M_PI * freq2) * x);
        data[i] = ScalarType(values[i], grads[i]);
    }
    VectorND<ComplexType> answer{};
    /* Make answer */ {
        FFT<ValueType> fft(N, PlanFlag::Estimate);
        fft.getRSpace() = values;
        fft.transform();
        VectorND<ComplexPlainScalar> k_values = fft.getKSpace();
        fft.getRSpace() = grads;
        fft.transform();
        VectorND<ComplexPlainScalar> k_grads = fft.getKSpace();
        if (k_values.getLength() != k_grads.getLength()) [[unlikely]]
            exit(EXIT_FAILURE);

        answer.resize(k_values.getLength());
        for (size_t i = 0; i < answer.getLength(); ++i)
            answer[i] = ComplexType(k_values[i], k_grads[i]);
    }
    FFT<ScalarType> fft(N, PlanFlag::Estimate);
    /* Test transform */ {
        fft.getRSpace() = data;
        fft.transform();
        if (!vectorNear(fft.getKSpace(), answer, 1E-15))
            exit(EXIT_FAILURE);
    }
    /* Test invTransform */ {
        constexpr double precision = 1E-13;
        fft.invTransform();
        for (size_t i = 0; i < data.getLength(); ++i) {
            const bool isNear = scalarNear(ScalarType(data[i]), ScalarType(fft.getRSpace()[i]), precision);
            const bool isSmall = abs(data[i].getValue()) < ValueType(precision) && abs(fft.getRSpace()[i].getValue()) < ValueType(precision);
            if(!isNear && !isSmall)
                exit(EXIT_FAILURE);
        }
    }
}

int main() {
    /* 1D real */ {
        const size_t N = 100;
        const double t_max = 2;
        constexpr double freq1 = 3;
        constexpr double freq2 = 4;
        
        VectorND<RealType> data(N);
        {
            const VectorND<RealType> v_x = VectorND<RealType>::linspace(RealType(0), RealType(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(RealType(2 * M_PI * freq1) * x) + sin(RealType(2 * M_PI * freq2) * x) * 2;
            }
        }
        FFT<RealType> fft(data, PlanFlag::Measure);
        const VectorND<RealType> intense = toNormVector(fft.getKSpace());

        /* Parseval theorem */ {
            const RealType power = square(data).sum();
            const RealType power_fft = square(intense).sum() / RealType(intense.getLength() - 1);
            if (!scalarNear(power, power_fft, 1E-15))
                return 1;
        }
        /* Test freq */ {
            const double deltaFreq = double(fft.getKSpaceDelta(RealType(t_max / N))) / (2 * M_PI);
            const RealType freq1_power = intense[freq1 / deltaFreq];
            const RealType freq2_power = intense[freq2 / deltaFreq];
            if (!scalarNear(freq2_power / freq1_power, RealType(2), 1E-14))
                return 1;
        }
        /* Test inv */ {
            constexpr double precision = 1E-13;
            fft.invTransform();
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
        
        VectorND<ComplexType> data(N);
        {
            const VectorND<RealType> v_x = VectorND<RealType>::linspace(RealType(0), RealType(t_max), N + 1);
            for (size_t i = 0; i < N; ++i) {
                const auto& x = v_x[i];
                data[i] = sin(RealType(2 * M_PI * freq1) * x) + sin(RealType(2 * M_PI * freq2) * x) * 2;
            }
        }
        FFT<ComplexType> fft(data, PlanFlag::Estimate);
        VectorND<ComplexType> trans(N);
        for (size_t i = 0; i < N; ++i) {
            ComplexType temp(0);
            for (size_t j = 0; j < N; ++j) {
                RealType phase = 2 * M_PI * j * i / N;
                temp += data[j] * ComplexType(cos(phase), -sin(phase));
            }
            trans[i] = temp;
        }
        fft.invTransform();
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

        DenseMatrix<RealType, MatrixOption::Row | MatrixOption::Vector> data(N1, N2);
        for (size_t i = 0; i < N1; ++i)
            for (size_t j = 0; j < N2; ++j)
                data(i, j) = RealType(std::sin(2 * M_PI * freq1 * i * deltaX) + 2 * std::cos(2 * M_PI * freq2 * j * deltaY));

        FFT<RealType, 2> fft({N1, N2}, PlanFlag::Estimate);
        fft.transform(data);
        /* Test freq */ {
            const double deltaFreq1 = double(fft.getKSpaceDelta(deltaX, 0)) / (2 * M_PI);
            const double deltaFreq2 = double(fft.getKSpaceDelta(deltaY, 1)) / (2 * M_PI);
            const VectorND<RealType> intense = toNormVector(fft.getKSpace().flatten());
            const RealType freq1_power = intense[size_t(freq1 / deltaFreq1) * fft.getKSpaceSize()[1]];
            const RealType freq1_power_conj = intense[(N1 - size_t(freq1 / deltaFreq1)) * fft.getKSpaceSize()[1]];
            const RealType freq2_power = intense[freq2 / deltaFreq2];
            if (!scalarNear(freq1_power, freq1_power_conj, 1E-15))
                return 1;
            if (!scalarNear(freq2_power / freq1_power, RealType(2), 1E-14))
                return 1;
        }
        /* Test inv */ {
            constexpr double precision = 1E-11;
            fft.invTransform();
            for (size_t i = 0; i < data.getRow(); ++i) {
                for (size_t j = 0; j < data.getCol(); ++j) {
                    const bool isNear = scalarNear(data(i, j), fft.getRSpace()(i, j), precision);
                    const bool isSmall = abs(data(i, j)) < RealType(precision) && abs(fft.getRSpace()(i, j)) < RealType(precision);
                    if(!isNear && !isSmall)
                        return 1;
                }
            }
        }
    }
    test_differentiable();
    return 0;
}
