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
        template<Scalar T, ExprType Type>
        __global__ void __launch_bounds__(1, 1) ElementaryFunction_calcKernel(
                Physica::PlainStruct<device_obj<TraceSegment<T, 1>>> segment_,
                Physica::PlainStruct<const device_obj<Diff<T, DiffMode::Reverse, 1>>> s_) {
            using SegmentType = device_obj<TraceSegment<T, 1>>;
            using DiffRecord = typename SegmentType::DiffRecord;
            auto& segment = segment_.getDerived();
            segment.getRecords()[0] = DiffRecord{0, Type};
            const auto& s = s_.getDerived();
            segment.getOperands()[0] = s;
            if constexpr (Type == ExprType::Ln)
                segment.getValues()[0] = ln(s.getValue());
            else
                static_assert(Type == ExprType::Ln, "[Error]: Not implemented");
            segment.getGrads()[0] = 0;
        }
    }

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> ln(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s) {
        using TracerType = device_obj<DiffTracer<T, 1>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExprType::Ln);
        Internal::ElementaryFunction_calcKernel<T, ExprType::Ln>
                <<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(segment), asStruct(s));
        check(cudaGetLastError());
        return segment[0];
    }
}
