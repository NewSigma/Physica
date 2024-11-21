/*
 * Copyright 2024 Weibo He.
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
        template<Matrix T, Vector V>
        __global__ void DiffMatrixProduct_gemvKernel(
                Physica::PlainStruct<const T> m_,
                Physica::PlainStruct<const V> v_,
                Physica::PlainStruct<T> dots_,
                Physica::PlainStruct<V> result_) {
            using DiffRecord = V::DiffRecord;
            using ScalarType = V::ScalarType;
            using ValueType = ScalarType::ValueType;
            extern __shared__ ValueType buffer[];
            const T& m = m_.getDerived();
            const V& v = v_.getDerived();
            T& dots = dots_.getDerived();
            const int col = m.getCol();

            const unsigned int r = blockIdx.x;
            const int index = threadIdx.x;
            ValueType threadSum = 0;
            for (int i = index; i < col; i += blockDim.x) {
                const size_t offset = dots.calcOffset(r, i);
                dots.getRecord(r, i) = DiffRecord{offset * 2, ExprType::Mul};
                dots.getOperands()[2 * offset] = m.calc(r, i);
                dots.getOperands()[2 * offset + 1] = v.calc(i);
                ValueType value = m.getValue(r, i) * v.getValue(i);
                threadSum += value;
                dots.getValue(r, i) = value;
                dots.getGrad(r, i) = 0;
            }
            buffer[index] = threadSum;

            int i0 = blockDim.x;
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
            V& result = result_.getDerived();
            result.getRecord(r) = DiffRecord{r * 2, ExprType::Sum};
            result.getOperands()[r * 2] = dots.calc(r, 0);
            result.getOperands()[r * 2 + 1] = dots.calc(r, col - 1);
            result.getValue(r) = threadSum;
            result.getGrad(r) = 0;
        }
    }

    template<Scalar T, int Option, int Order>
    auto operator*(
            const device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>& m,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v) {
        using M = device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>;
        using V = device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>;
        assert(m.getRow() > 0 && "[Error]: This is a empty matrix");
        assert(m.getCol() == v.getLength() && "[Error]: Dims do not match");
        M dots(m.getRow(), m.getCol(), ExprType::Mul);
        V result(m.getRow(), ExprType::Sum);
        const size_t numThread = std::min(m.getCol(), V::MaxThreadPerBlock);
        Internal::DiffMatrixProduct_gemvKernel<M, V>
                <<<m.getRow(), numThread, numThread * sizeof(T), CUDAContext::getInstance()>>>(asStruct(m), asStruct(v), asStruct(dots), asStruct(result));
        check(cudaGetLastError());
        return result;
    }
}
