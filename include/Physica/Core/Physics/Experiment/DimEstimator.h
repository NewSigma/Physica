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

#include <algorithm>
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/LinearFit.h"
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    /**
     * References:
     * [1] Neural Process. Lett. 14, 27–34 (2001). https://doi.org/10.1023/A:1011326007550
     * [2] Physica D 9, 189-208 (1983); https://doi.org/10.1016/0167-2789(83)90298-1
     */
    class DimEstimator {
        using ScalarType = float64;
        using DataMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector>;
        VectorND<ScalarType> intrinsicDim;
        VectorND<ScalarType> correlateDim;
    public:
        template<Vector T, RNG R>
        DimEstimator(size_t sampleNum, const Array<size_t>& intrinsicDim_, const T& radius);
        /* Operations */
        template<Vector T, Matrix U>
        ScalarType intrinDim(const U& data, const T& radius) const;
        template<Vector T, Matrix U>
        static ScalarType corrDimen(const U& data, const T& radius);
        template<Vector T>
        static DataMatrix toHighDimForm(const T& data, size_t step, size_t dim);
    private:
        template<Vector T, Matrix U>
        static inline T corrIntegral(const U& data, const T& radius);
    };

    template<Vector T, RNG R>
    DimEstimator::DimEstimator(size_t sampleNum,
                               const Array<size_t>& intrinsicDim_,
                               const T& radius)
            : intrinsicDim(intrinsicDim_.getLength()) {
        const size_t length = intrinsicDim.getLength();
        correlateDim.resize(length);
        const size_t maxDim = *std::max_element(intrinsicDim_.cbegin(), intrinsicDim_.cend());
        const auto data = DataMatrix::template random_uniform<R>(sampleNum, maxDim);
        for (size_t i = 0; i < length; ++i) {
            intrinsicDim[i] = ScalarType(intrinsicDim_[i]);
            correlateDim[i] = corrDimen(data.leftCols(intrinsicDim_[i]), radius);
        }
    }

    template<Vector T, Matrix U>
    DimEstimator::ScalarType
    DimEstimator::intrinDim(const U& data, const T& radius) const {
        return lagrange(intrinsicDim, correlateDim, corrDimen(data, radius));
    }
    /**
     * \param data
     * Each row represents a piece of data
     */
    template<Vector T, Matrix U>
    DimEstimator::ScalarType
    DimEstimator::corrDimen(const U& data, const T& radius) {
        const T logCorrIntegral = ln(corrIntegral(data, radius));
        const T logR = ln(radius);
        return LinearFit<ScalarType>::fit(logR, logCorrIntegral).first;
    }
    /**
     * Helper function for distinguishing chaos and random noise, refer to [2]
     */
    template<Vector T>
    DimEstimator::DataMatrix DimEstimator::toHighDimForm(
            const T& data,
            size_t step,
            size_t dim) {
        assert(step > 0);
        assert(dim > 1);
        const size_t newNumData = (data.getLength() - (dim - 1) * step) / step;
        DataMatrix result(newNumData, dim);
        for (size_t i = 0; i < newNumData; ++i) {
            auto row = result.row(i);
            for (size_t j = 0; j < dim; ++j)
                row[j] = data[i + step * j];
        }
        return result;
    }

    template<Vector T, Matrix U>
    inline T DimEstimator::corrIntegral(const U& data, const T& radius) {
        const size_t numData = data.getRow();
        const T squaredR = square(radius);
        T count(radius.getLength(), 0);
        for (size_t i = 0; i < numData - 1; ++i) {
            auto data1 = data.row(i);
            for (size_t j = i + 1; j < numData; ++j) {
                auto data2 = data.row(j);
                const ScalarType r2 = (data1 - data2).squaredNorm();
                for (size_t k = 0; k < count.getLength(); ++k)
                    count[k] += ScalarType(r2 <= squaredR[k]);
            }
        }
        const ScalarType factor = ScalarType(2) / ScalarType(numData * (numData - 1));
        return factor * count;
    }
}
