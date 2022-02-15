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

#include <algorithm>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/LValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorExpression.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/LinearFit.h"
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    /**
     * References:
     * [1] F. Camastra and A. Vinciarelli, Neural Process. Lett. 14, 27 (2001).
     * [2] Grassberger, P. and Procaccia, I.: Measuring the strangeness of strange attractors, Physica, D9 (1983) 189^208.
     */
    class DimEstimator {
        using ScalarType = Scalar<Double, false>;
        using DataMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector>;
        Vector<ScalarType> intrinsicDim;
        Vector<ScalarType> correlateDim;
    public:
        template<class VectorType>
        DimEstimator(size_t sampleNum, const Utils::Array<size_t>& intrinsicDim_, const LValueVector<VectorType>& radius);
        /* Operations */
        template<class MatrixType, class VectorType>
        ScalarType intrinDim(const LValueMatrix<MatrixType>& data, const LValueVector<VectorType>& radius);
        template<class MatrixType, class VectorType>
        static ScalarType corrDimen(const LValueMatrix<MatrixType>& data, const LValueVector<VectorType>& radius);
        template<class VectorType>
        static DataMatrix toHighDimForm(const LValueVector<VectorType>& data, size_t step, size_t dim);
    private:
        template<class MatrixType, class VectorType>
        static inline VectorType corrIntegral(const LValueMatrix<MatrixType>& data, const LValueVector<VectorType>& radius);
    };

    template<class VectorType>
    DimEstimator::DimEstimator(size_t sampleNum,
                               const Utils::Array<size_t>& intrinsicDim_,
                               const LValueVector<VectorType>& radius)
            : intrinsicDim(intrinsicDim_.getLength()) {
        const size_t length = intrinsicDim.getLength();
        correlateDim.resize(length);
        const size_t maxDim = *std::max_element(intrinsicDim_.cbegin(), intrinsicDim_.cend());
        const auto data = DataMatrix::randomMatrix(sampleNum, maxDim);
        for (size_t i = 0; i < length; ++i) {
            intrinsicDim[i] = ScalarType(intrinsicDim_[i]);
            correlateDim[i] = corrDimen(data.leftCols(intrinsicDim_[i]), radius);
        }
    }

    template<class MatrixType, class VectorType>
    typename DimEstimator::ScalarType
    DimEstimator::intrinDim(const LValueMatrix<MatrixType>& data, const LValueVector<VectorType>& radius) {
        return lagrange(intrinDim, correlateDim, corrDimen(data, radius));
    }
    /**
     * \param data
     * Each row represents a piece of data
     */
    template<class MatrixType, class VectorType>
    typename DimEstimator::ScalarType
    DimEstimator::corrDimen(const LValueMatrix<MatrixType>& data, const LValueVector<VectorType>& radius) {
        const VectorType logCorrIntegral = ln(corrIntegral(data, radius));
        const VectorType logR = ln(radius);
        return linearFit(logR, logCorrIntegral).first;
    }
    /**
     * Helper function for distinguishing chaos and random noise, refer to [2]
     */
    template<class VectorType>
    typename DimEstimator::DataMatrix DimEstimator::toHighDimForm(
            const LValueVector<VectorType>& data,
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

    template<class MatrixType, class VectorType>
    inline VectorType DimEstimator::corrIntegral(const LValueMatrix<MatrixType>& data, const LValueVector<VectorType>& radius) {
        const size_t numData = data.getRow();
        const VectorType squaredR = square(radius);
        VectorType count(radius.getLength(), 0);
        for (size_t i = 0; i < numData - 1; ++i) {
            auto data1 = data.row(i);
            for (size_t j = i + 1; j < numData; ++j) {
                auto data2 = data.row(j);
                count += (data1.asVector() - data2).squaredNorm() <= squaredR;
            }
        }
        const ScalarType factor = ScalarType::Two() / ScalarType(numData * (numData - 1));
        return factor * count;
    }
}
