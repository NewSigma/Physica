/*
 * Copyright 2024 WeiBo He.
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

namespace Physica::Core {
    namespace Internal {
        template<class MatrixType, class VectorType>
        __global__ void DiffMatrixProduct_gemvKernel(
                Physica::PlainStruct<const MatrixType> m_,
                Physica::PlainStruct<const VectorType> v_,
                Physica::PlainStruct<MatrixType> dots_,
                Physica::PlainStruct<VectorType> result_) {
            using DiffRecord = typename VectorType::DiffRecord;
            using ScalarType = typename VectorType::ScalarType;
            using PlainScalar = typename ScalarType::PlainScalar;
            extern __shared__ PlainScalar buffer[];
            const MatrixType& m = m_.getDerived();
            const VectorType& v = v_.getDerived();
            MatrixType& dots = dots_.getDerived();
            const int column = m.getColumn();

            const unsigned int r = blockIdx.x;
            const int index = threadIdx.x;
            PlainScalar threadSum = 0;
            for (int i = index; i < column; i += blockDim.x) {
                const size_t offset = dots.calcOffset(r, i);
                dots.getRecord(r, i) = DiffRecord{offset * 2, ExpressionType::Mul};
                dots.getOperands()[offset] = m.calc(r, i);
                dots.getOperands()[offset + 1] = v.calc(i);
                PlainScalar value = m.getValue(r, i) * v.getValue(i);
                threadSum += value;
                dots.getValue(r, i) = value;
                dots.getGrad(r, i) = 0;
            }
            buffer[index] = threadSum;

            int i0 = column;
            int i = (i0 + 1) / 2;
            while (index + i < i0) {
                __syncthreads();
                threadSum += buffer[index + i];
                buffer[index] = threadSum;
                i0 = i;
                i = (i0 + 1) / 2;
            }

            if (index != 0)
                return;
            VectorType& result = result_.getDerived();
            result.getRecord(r) = DiffRecord{r * 2, ExpressionType::Sum};
            result.getOperands()[r * 2] = dots.calc(r, 0);
            result.getOperands()[r * 2 + 1] = dots.calc(r, column - 1);
            result.getValue(r) = threadSum;
            result.getGrad(r) = 0;
        }
    }

    template<class PlainScalar, int Option>
    auto operator*(
            const device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>>& m,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v) {
        using MatrixType = device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>>;
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>;
        assert(m.getRow() > 0 && "[Error]: This is a empty matrix");
        assert(m.getColumn() == v.getLength() && "[Error]: Dims do not match");
        MatrixType dots(m.getRow(), m.getColumn(), ExpressionType::Mul);
        VectorType result(m.getRow(), ExpressionType::Sum);
        const size_t numThread = m.getColumn() > VectorType::MaxThreadPerBlock ? VectorType::MaxThreadPerBlock : m.getColumn();
        Internal::DiffMatrixProduct_gemvKernel<MatrixType, VectorType>
                <<<m.getRow(), numThread, numThread * sizeof(PlainScalar), StreamPool::getStream()>>>(asStruct(m), asStruct(v), asStruct(dots), asStruct(result));
        return result;
    }
}
