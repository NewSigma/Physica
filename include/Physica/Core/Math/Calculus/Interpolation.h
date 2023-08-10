/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica::Core {
    namespace Internal {
        /**
         * Critical to performance of \class Ewald, refactor with care
         */
        template<class ScalarType>
        inline ScalarType quadraticInterpolate(ScalarType x1, ScalarType x2, ScalarType x3, ScalarType y1, ScalarType y2, ScalarType y3, ScalarType x) {
            const ScalarType x_x1 = x - x1;
            const ScalarType x_x2 = x - x2;
            const ScalarType x_x3 = x - x3;
            const ScalarType x1_x2 = x1 - x2;
            const ScalarType x2_x3 = x2 - x3;
            const ScalarType x3_x1 = x3 - x1;
            return -((x_x2 * x_x3) * (x2_x3 * y1) + (x_x1 * x_x3) * (x3_x1 * y2) + (x_x1 * x_x2) * (x1_x2 * y3)) / (x1_x2 * x2_x3 * x3_x1);
        }

        template<class ScalarType>
        inline ScalarType quadraticInterpolate_diff(ScalarType x1, ScalarType x2, ScalarType x3, ScalarType y1, ScalarType y2, ScalarType y3, ScalarType x) {
            const ScalarType xx = x * 2.0;
            const ScalarType x1_xx = x1 - xx;
            const ScalarType x2_xx = x2 - xx;
            const ScalarType x3_xx = x3 - xx;
            const ScalarType y1_y2 = y1 - y2;
            const ScalarType y2_y3 = y2 - y3;
            const ScalarType y3_y1 = y3 - y1;
            const ScalarType x1_x2 = x1 - x2;
            const ScalarType x2_x3 = x2 - x3;
            const ScalarType x3_x1 = x3 - x1;
            return (x3 * x3_xx * y1_y2 + x1 * x1_xx * y2_y3 + x2 * x2_xx * y3_y1) / (x1_x2 * x2_x3 * x3_x1);
        }
        /**
         * \param factor
         * Equals to ((x1 - x2) * (x2 - x3) * (x3 - x1))
         */
        template<class ScalarType>
        inline ScalarType quadraticInterpolate_diff1(ScalarType factor, ScalarType step, ScalarType x2, ScalarType y1, ScalarType y2, ScalarType y3, ScalarType x) {
            const ScalarType repFactor = reciprocal(factor);
            const ScalarType y1_y2 = y1 - y2;
            const ScalarType y2_y3 = y2 - y3;
            const ScalarType y3_y1 = y3 - y1;
            return -(step * y3_y1 + ScalarType(2) * (x - x2) * (y1_y2 - y2_y3)) * repFactor;
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009:81-82
     */
    template<class VectorType>
    typename VectorType::ScalarType lagrange(
            const LValueVector<VectorType>& x0,
            const LValueVector<VectorType>& y0,
            typename VectorType::ScalarType x) {
        using ScalarType = typename VectorType::ScalarType;
        assert(x0.getLength() == y0.getLength());

        const size_t length = x0.getLength();
        VectorType c = y0;
        VectorType d = y0;
        size_t ns = 0;
        /* Init ns */ {
            ScalarType delta = abs(x - x0[0]);
            for (size_t i = 1; i < length; ++i) {
                const ScalarType temp = abs(x - x0[i]);
                if (temp < delta) {
                    ns = i;
                    delta = temp;
                }
            }
        }
        ScalarType result = y0[ns--];
        for (size_t i = 1; i < length; ++i) {
            const size_t max_j = length - i;
            for (size_t j = 0; j < max_j; ++j) {
                const ScalarType delta1 = x0[j] - x;
                const ScalarType delta2 = x0[j + i] - x;
                const ScalarType factor = (c[j + 1] - d[j]) / (delta1 - delta2);
                c[j] = delta1 * factor;
                d[j] = delta2 * factor;
            }
            result += (2 * (ns + 1) < max_j) ? c[ns + 1] : d[ns--];
        }
        return result;
    }

    template<class VectorType>
    VectorType interpolate_fft(const LValueVector<VectorType>& data, size_t newDim) {
        using ScalarType = typename VectorType::ScalarType;
        using ComplexType = typename ScalarType::ComplexType;

        const size_t rSpaceSize = data.getLength();
        FFT<ScalarType, 1> fft(rSpaceSize, 1, PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();
        const size_t kSpaceSize = fft.getKSpaceSize();
        const auto& kSpace = fft.getKSpace();

        auto result = VectorType(newDim);
        for (size_t i = 0; i < newDim; ++i) {
            auto elem = ComplexType(0);
            for (size_t k = 0; k < rSpaceSize; ++k) {
                const ScalarType phase = 2 * M_PI * (k * i) / newDim;
                ScalarType s, c;
                sincos(phase, s, c);
                const auto factor = ComplexType(c, s);
                if (k >= kSpaceSize)
                    elem += kSpace[rSpaceSize - k].conjugate() * factor;
                else
                    elem += kSpace[k] * factor;
            }

            if constexpr (ScalarType::isComplex)
                result[i] = elem;
            else
                result[i] = elem.getReal();
        }
        result *= reciprocal(ScalarType(rSpaceSize));
        return result;
    }

    template<class VectorType>
    typename VectorType::ScalarType interpolate_fft(
            const LValueVector<VectorType>& data,
            typename VectorType::ScalarType x,
            typename VectorType::ScalarType period) {
        using ScalarType = typename VectorType::ScalarType;
        using ComplexType = typename ScalarType::ComplexType;

        const size_t rSpaceSize = data.getLength();
        FFT<ScalarType, 1> fft(rSpaceSize, 1, PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();
        const size_t kSpaceSize = fft.getKSpaceSize();
        const auto& kSpace = fft.getKSpace();

        const ScalarType relative_x = x / period;
        auto elem = ComplexType(0);
        for (size_t k = 0; k < rSpaceSize; ++k) {
            const ScalarType phase = ScalarType(2 * M_PI * k) * relative_x;
            ScalarType s, c;
            sincos(phase, s, c);
            const auto factor = ComplexType(c, s);
            if (k >= kSpaceSize)
                elem += kSpace[rSpaceSize - k].conjugate() * factor;
            else
                elem += kSpace[k] * factor;
        }

        ScalarType result;
        if constexpr (ScalarType::isComplex)
            result = elem;
        else
            result = elem.getReal();
        result *= reciprocal(ScalarType(rSpaceSize));
        return result;
    }

    template<class GridType>
    GridType interpolate_fft(const LValueGrid<GridType>& data, Utils::Array<size_t, 3> newDim) {
        using ScalarType = typename GridType::ScalarType;
        using ComplexType = typename ScalarType::ComplexType;
        using Index3D = Utils::Array<size_t, 3>;
        constexpr int Dim = 3;

        FFT<ScalarType, 3> fft(data.getDim(), 1, PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();

        auto result = GridType(newDim);
        result.forIndexInGrid([newDim, &fft, &result](Index3D rIndex) {
            const auto& kSpace = fft.getKSpace();
            auto elem = ComplexType(0);
            kSpace.forIndexInGrid([newDim, rIndex, &fft, &kSpace, &elem](Index3D kIndex) {
                const Index3D rSpaceSize = fft.getRSpaceSize();
                const Index3D kSpaceSize = fft.getKSpaceSize();
                ScalarType phase = 0, s, c;
                for (int dim = 0; dim < Dim; ++dim)
                    phase += ScalarType(kIndex[dim] * rIndex[dim]) / newDim[dim];
                phase *= ScalarType(2 * M_PI);
                sincos(phase, s, c);
                const auto factor = ComplexType(c, s);
                if (kIndex[2] >= kSpaceSize[2])
                    elem += kSpace[rSpaceSize[2] - kIndex[2]].conjugate() * factor;
                else
                    elem += kSpace[kIndex[2]] * factor;
            });

            if constexpr (ScalarType::isComplex)
                result(rIndex) = elem;
            else
                result(rIndex) = elem.getReal();
        });
        result *= reciprocal(reciprocal(data.getSize()));
        return result;
    }

    template<class GridType>
    typename GridType::ScalarType interpolate_fft(
            const LValueGrid<GridType>& data,
            Vector<typename GridType::ScalarType, 3> r,
            Vector<typename GridType::ScalarType, 3> period) {
        using ScalarType = typename GridType::ScalarType;
        using ComplexType = typename ScalarType::ComplexType;
        using Index3D = Utils::Array<size_t, 3>;
        constexpr int Dim = 3;
        
        FFT<ScalarType, 3> fft(data.getDim(), 1, PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();
        const auto& kSpace = fft.getKSpace();

        const Vector<ScalarType, 3> relative_r = r / period;
        auto elem = ComplexType(0);
        kSpace.forIndexInGrid([relative_r, &fft, &kSpace, &elem](Index3D kIndex) {
            const Index3D rSpaceSize = fft.getRSpaceSize();
            const Index3D kSpaceSize = fft.getKSpaceSize();
            ScalarType phase = 0, s, c;
            for (int dim = 0; dim < Dim; ++dim)
                phase += ScalarType(kIndex[dim]) * relative_r[dim];
            phase *= ScalarType(2 * M_PI);
            sincos(phase, s, c);
            const auto factor = ComplexType(c, s);
            if (kIndex[2] >= kSpaceSize[2])
                elem += kSpace[rSpaceSize[2] - kIndex[2]].conjugate() * factor;
            else
                elem += kSpace[kIndex[2]] * factor;
        });

        ScalarType result;
        if constexpr (ScalarType::isComplex)
            result = elem;
        else
            result = elem.getReal();
        result *= reciprocal(reciprocal(data.getSize()));
        return result;
    }
}
