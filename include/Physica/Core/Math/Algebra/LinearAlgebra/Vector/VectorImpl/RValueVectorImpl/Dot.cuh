/*
 * Copyright 2024-2026 Weibo He.
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

#include "../RValueVector.cuh"

namespace Physica {
    template<Vector V1, Vector V2>
    class device_obj<Dot<V1, V2>> {
        using This = device_obj<Dot<V1, V2>>;
        using T1 = V1::ScalarType;
        using T2 = V2::ScalarType;
        using T = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        using Tr = T::RealType;
    private:
        const V1& v1;
        const V2& v2;
    public:
        __host__ __device__ device_obj(const V1& v1_, const V2& v2_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ T calc() const noexcept;
        [[nodiscard]] __device__ T calc(instanceof_x<ThreadBlock> auto block) const noexcept;
    };

    template<Vector V1, Vector V2>
    __host__ __device__ device_obj<Dot<V1, V2>>::device_obj(const V1& v1_, const V2& v2_) : v1(v1_), v2(v2_) {
        assert(v1.getLength() == v2.getLength());
    }

    template<Vector V1, Vector V2>
    __host__ __device__ auto device_obj<Dot<V1, V2>>::calc() const noexcept -> T {
        if (IsHost()) {
            auto kernel = [v1_ = asStruct(v1), v2_ = asStruct(v2)] __device__() mutable {
                const auto& v1 = v1_.getDerived();
                const auto& v2 = v2_.getDerived();

                ThreadBlock<CUDADevAttr::DefaultThreadsPerBlock> block{};
                T local = 0;
                for (size_t i = block.tid(); i < v1.getLength(); i += block.getNumThread()) {
                    if constexpr (std::same_as<T1, T2>)
                        local = fma(v1.calc(i), v2.calc(i), local);
                    else
                        local += v1.calc(i) * v2.calc(i);
                }
                return block.sum(local);
            };
            return CUDAExecutor::launch(kernel, KernelConfig(1, CUDADevAttr::DefaultThreadsPerBlock));
        }

        if constexpr (IsDevice())
            return calc(ThreadBlock<1>{});
    }

    template<Vector V1, Vector V2>
    __device__ auto device_obj<Dot<V1, V2>>::calc(instanceof_x<ThreadBlock> auto block) const noexcept -> T {
        T dot = 0;
        const size_t length = v1.getLength();
        for (size_t i = block.tid(); i < length; i += block.getNumThread()) {
            if constexpr (std::same_as<T1, T2>)
                dot = fma(v1.calc(i), v2.calc(i), dot);
            else
                dot += v1.calc(i) * v2.calc(i);
        }
        return block.sum(dot);
    }

    template<Vector V1, Vector V2>
    __host__ __device__ auto dot(const V1& v1, const V2& v2) requires(DeviceObj<V1> && DeviceObj<V2>) {
        if constexpr (!canonicalized(v1, v2))
            return device_obj<Dot<V2, V1>>(v2, v1);
        else
            return device_obj<Dot<V1, V2>>(v1, v2);
    }
}
