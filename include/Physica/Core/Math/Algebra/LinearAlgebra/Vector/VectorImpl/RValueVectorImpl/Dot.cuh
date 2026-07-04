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
    template<Vector LHS, Vector RHS>
    class device_obj<Dot<LHS, RHS>> {
        using This = device_obj<Dot<LHS, RHS>>;
        using T1 = std::remove_cvref_t<LHS>::ScalarType;
        using T2 = std::remove_cvref_t<RHS>::ScalarType;
        using T = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        using Tr = T::RealType;
    protected:
        using Ref1 = add_device_obj<LHS>::type;
        using Ref2 = add_device_obj<RHS>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_cvref_t<LHS>>> lhs;
        PlainStruct<add_device_obj_t<std::remove_cvref_t<RHS>>> rhs;
    public:
        __host__ __device__ device_obj(Ref1 lhs_, Ref2 rhs_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ T calc() const noexcept;
        [[nodiscard]] __device__ T calc(instanceof_x<ThreadBlock> auto block) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ constexpr auto&& getRHS(this auto&&) noexcept;
    };

    template<Vector LHS, Vector RHS>
    __host__ __device__ device_obj<Dot<LHS, RHS>>::device_obj(Ref1 lhs_, Ref2 rhs_) : lhs(asStruct(lhs_)), rhs(asStruct(rhs_)) {
        assert(getLHS().getLength() == getRHS().getLength());
    }

    template<Vector LHS, Vector RHS>
    __host__ __device__ auto device_obj<Dot<LHS, RHS>>::calc() const noexcept -> T {
        if (IsHost()) {
            auto kernel = [lhs_ = lhs, rhs_ = rhs] __device__() mutable {
                const auto& lhs = lhs_.getDerived();
                const auto& rhs = rhs_.getDerived();

                ThreadBlock<CUDADevAttr::DefaultThreadsPerBlock> block{};
                T local = 0;
                for (size_t i = block.tid(); i < lhs.getLength(); i += block.getNumThread()) {
                    if constexpr (std::same_as<T1, T2>)
                        local = fma(lhs.calc(i), rhs.calc(i), local);
                    else
                        local += lhs.calc(i) * rhs.calc(i);
                }
                return block.sum(local);
            };
            return CUDAExecutor::launch(kernel, KernelConfig(1, CUDADevAttr::DefaultThreadsPerBlock));
        }

        if constexpr (IsDevice())
            return calc(ThreadBlock<1>{});
    }

    template<Vector LHS, Vector RHS>
    __device__ auto device_obj<Dot<LHS, RHS>>::calc(instanceof_x<ThreadBlock> auto block) const noexcept -> T {
        T dot = 0;
        const size_t length = getLHS().getLength();
        for (size_t i = block.tid(); i < length; i += block.getNumThread()) {
            if constexpr (std::same_as<T1, T2>)
                dot = fma(getLHS().calc(i), getRHS().calc(i), dot);
            else
                dot += getLHS().calc(i) * getRHS().calc(i);
        }
        return block.sum(dot);
    }

    template<Vector LHS, Vector RHS>
    __host__ __device__ constexpr auto&& device_obj<Dot<LHS, RHS>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref1>(self.lhs.getDerived());
    }

    template<Vector LHS, Vector RHS>
    __host__ __device__ constexpr auto&& device_obj<Dot<LHS, RHS>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref2>(self.rhs.getDerived());
    }

    template<Vector LHS, Vector RHS>
    __host__ __device__ auto dot(LHS&& lhs, RHS&& rhs) noexcept requires(DeviceObj<LHS> && DeviceObj<RHS>) {
        if constexpr (!canonicalized(lhs, rhs))
            return device_obj<Dot<remove_device_obj_t<RHS&&>, remove_device_obj_t<LHS&&>>>(std::forward<RHS>(rhs), std::forward<LHS>(lhs));
        else
            return device_obj<Dot<remove_device_obj_t<LHS&&>, remove_device_obj_t<RHS&&>>>(std::forward<LHS>(lhs), std::forward<RHS>(rhs));
    }
}
