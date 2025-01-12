/*
 * Copyright 2022-2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica::Core {
    namespace Internal {
        /**
         * Critical to performance of \class Ewald, refactor with care
         */
        template<Scalar T>
        __host__ __device__ inline T quadraticInterpolate(
                T x1, T x2, T x3, T y1, T y2, T y3, T x) {
            const T x_x1 = x - x1;
            const T x_x2 = x - x2;
            const T x_x3 = x - x3;
            const T x1_x2 = x1 - x2;
            const T x2_x3 = x2 - x3;
            const T x3_x1 = x3 - x1;
            const T factor = x1_x2 * x2_x3 * x3_x1;
            return -((x_x2 * x_x3) * (x2_x3 * y1) + (x_x1 * x_x3) * (x3_x1 * y2) + (x_x1 * x_x2) * (x1_x2 * y3)) / factor;
        }

        template<Scalar T>
        __host__ __device__ inline T quadraticInterpolate_diff(
                T x1, T x2, T x3, T y1, T y2, T y3, T x) {
            const T xx = x * 2.0;
            const T x1_xx = x1 - xx;
            const T x2_xx = x2 - xx;
            const T x3_xx = x3 - xx;
            const T y1_y2 = y1 - y2;
            const T y2_y3 = y2 - y3;
            const T y3_y1 = y3 - y1;
            const T x1_x2 = x1 - x2;
            const T x2_x3 = x2 - x3;
            const T x3_x1 = x3 - x1;
            const T factor = x1_x2 * x2_x3 * x3_x1;
            return -(x3 * x3_xx * y1_y2 + x1 * x1_xx * y2_y3 + x2 * x2_xx * y3_y1) / factor;
        }
        /**
         * \param factor
         * Equals to reciprocal((x1 - x2) * (x2 - x3) * (x3 - x1))
         */
        template<Scalar T>
        __host__ __device__ inline T quadraticInterpolate_diff1(
                T factor, T step, T x2, T y1, T y2, T y3, T x) {
            const T y1_y2 = y1 - y2;
            const T y2_y3 = y2 - y3;
            const T y3_y1 = y3 - y1;
            return (step * y3_y1 + T(2) * (x - x2) * (y1_y2 - y2_y3)) * factor;
        }
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:81-82
     */
    template<Vector V>
    V::ScalarType lagrange(const V& x0, const V& y0, typename V::ScalarType x) {
        using T = V::ScalarType;
        assert(x0.getLength() == y0.getLength());

        const size_t length = x0.getLength();
        V c = y0;
        V d = y0;
        size_t ns = 0;
        /* Init ns */ {
            T delta = abs(x - x0[0]);
            for (size_t i = 1; i < length; ++i) {
                const T temp = abs(x - x0[i]);
                if (temp < delta) {
                    ns = i;
                    delta = temp;
                }
            }
        }
        T result = y0[ns--];
        for (size_t i = 1; i < length; ++i) {
            const size_t max_j = length - i;
            for (size_t j = 0; j < max_j; ++j) {
                const T delta1 = x0[j] - x;
                const T delta2 = x0[j + i] - x;
                const T factor = (c[j + 1] - d[j]) / (delta1 - delta2);
                c[j] = delta1 * factor;
                d[j] = delta2 * factor;
            }
            result += (2 * (ns + 1) < max_j) ? c[ns + 1] : d[ns--];
        }
        return result;
    }

    template<Vector V>
    V interpolate_fft(const V& data, size_t newDim) {
        using T = V::ScalarType;
        using RealType = T::RealType;
        using ComplexType = T::ComplexType;
        constexpr bool isComplex = T::isComplex;

        const size_t rSpaceSize = data.getLength();
        FFT<T, 1> fft(rSpaceSize, PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();
        const size_t kSpaceSize = fft.getKSpaceSize();
        const auto& kSpace = fft.getKSpace();

        auto result = V(newDim);
        for (size_t i = 0; i < newDim; ++i) {
            auto elem = ComplexType(0);
            for (size_t k = 0; k < rSpaceSize; ++k) {
                const RealType phase = 2 * M_PI * (k * i) / newDim;
                const auto factor = ComplexType::fromPhase(phase);
                if constexpr (!isComplex) {
                    if (k >= kSpaceSize)
                        elem += kSpace[rSpaceSize - k].conjugate() * factor;
                    else
                        elem += kSpace[k] * factor;
                }
                else
                    elem += kSpace[k] * factor;
            }

            if constexpr (isComplex)
                result[i] = elem;
            else
                result[i] = elem.real();
        }
        result *= reciprocal(RealType(rSpaceSize));
        return result;
    }

    template<Vector V>
    V::ScalarType interpolate_fft(const V& data, typename V::ScalarType::RealType x, typename V::ScalarType::RealType period) {
        using T = V::ScalarType;
        using RealType = T::RealType;
        using ComplexType = T::ComplexType;
        constexpr bool isComplex = T::isComplex;

        const size_t rSpaceSize = data.getLength();
        FFT<T, 1> fft(rSpaceSize, PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();
        const size_t kSpaceSize = fft.getKSpaceSize();
        const auto& kSpace = fft.getKSpace();

        const RealType relative_x = x / period;
        auto elem = ComplexType(0);
        for (size_t k = 0; k < rSpaceSize; ++k) {
            const RealType phase = RealType(2 * M_PI * k) * relative_x;
            const auto factor = ComplexType::fromPhase(phase);
            if constexpr (!isComplex) {
                if (k >= kSpaceSize)
                    elem += kSpace[rSpaceSize - k].conjugate() * factor;
                else
                    elem += kSpace[k] * factor;
            }
            else
                elem += kSpace[k] * factor;
        }

        T result;
        if constexpr (isComplex)
            result = elem;
        else
            result = elem.real();
        result *= reciprocal(RealType(rSpaceSize));
        return result;
    }

    template<LGrid G>
    RSpaceGrid<typename G::ScalarType> interpolate_fft(const G& data, typename GridBase::Index3D newDim) {
        using T = G::ScalarType;
        using ResultType = RSpaceGrid<T>;
        using RealType = T::RealType;
        using ComplexType = T::ComplexType;
        using Index3D = GridBase::Index3D;
        constexpr size_t Dim = 3;
        constexpr bool isComplex = T::isComplex;

        auto fft = FFT<T, 3>(data.getDim(), PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();

        auto result = ResultType(newDim);
        GridBase::forIndexInGrid(newDim, [newDim, &fft, &result](Index3D rIndex) {
            const auto& kSpace = fft.getKSpace();
            const Index3D rSpaceSize = fft.getRSpaceSize();
            auto elem = ComplexType(0);
            GridBase::forIndexInGrid(rSpaceSize, [newDim, rSpaceSize, rIndex, &kSpace, &elem](Index3D kIndex) {
                RealType phase(0);
                for (size_t dim = 0; dim < Dim; ++dim)
                    phase += RealType(kIndex[dim] * rIndex[dim]) / newDim[dim];
                phase *= RealType(2 * M_PI);
                const auto factor = ComplexType::fromPhase(phase);
                if constexpr (!isComplex) {
                    if (kIndex[2] >= kSpace.getDim()[2]) {
                        Index3D kIndex1;
                        kIndex1[0] = kIndex[0] == 0 ? kIndex[0] : (rSpaceSize[0] - kIndex[0]);
                        kIndex1[1] = kIndex[1] == 0 ? kIndex[1] : (rSpaceSize[1] - kIndex[1]);
                        kIndex1[2] = rSpaceSize[2] - kIndex[2];
                        elem += kSpace(kIndex1).conjugate() * factor;
                    }
                    else
                        elem += kSpace(kIndex) * factor;
                }
                else
                    elem += kSpace(kIndex) * factor;
            });

            if constexpr (isComplex)
                result(rIndex) = elem;
            else
                result(rIndex) = elem.real();
        });
        result *= reciprocal(RealType(data.getSize()));
        return result;
    }

    template<LGrid G>
    G::ScalarType interpolate_fft(
            const G& data,
            Vector3D<typename G::ScalarType::RealType> r,
            Vector3D<typename G::ScalarType::RealType> period) {
        using T = G::ScalarType;
        using RealType = T::RealType;
        using ComplexType = T::ComplexType;
        using Index3D = Array<size_t, 3>;
        constexpr size_t Dim = 3;
        constexpr bool isComplex = T::isComplex;
        
        auto fft = FFT<T, 3>(data.getDim(), PlanFlag::Estimate);
        fft.getRSpace() = data;
        fft.transform();
        const auto& kSpace = fft.getKSpace();

        const Vector3D<RealType> relative_r = divide(r, period);
        const Index3D rSpaceSize = fft.getRSpaceSize();
        auto elem = ComplexType(0);
        GridBase::forIndexInGrid(rSpaceSize, [relative_r, rSpaceSize, &kSpace, &elem](Index3D kIndex) {
            RealType phase(0);
            for (size_t dim = 0; dim < Dim; ++dim)
                phase += RealType(kIndex[dim]) * relative_r[dim];
            phase *= RealType(2 * M_PI);
            const auto factor = ComplexType::fromPhase(phase);
            if constexpr (!isComplex) {
                if (kIndex[2] >= kSpace.getDim()[2]) {
                    Index3D kIndex1;
                    kIndex1[0] = kIndex[0] == 0 ? kIndex[0] : (rSpaceSize[0] - kIndex[0]);
                    kIndex1[1] = kIndex[1] == 0 ? kIndex[1] : (rSpaceSize[1] - kIndex[1]);
                    kIndex1[2] = rSpaceSize[2] - kIndex[2];
                    elem += kSpace(kIndex1).conjugate() * factor;
                }
                else
                    elem += kSpace(kIndex) * factor;
            }
            else
                elem += kSpace(kIndex) * factor;
        });

        T result;
        if constexpr (isComplex)
            result = elem;
        else
            result = elem.real();
        result *= reciprocal(RealType(data.getSize()));
        return result;
    }
}
