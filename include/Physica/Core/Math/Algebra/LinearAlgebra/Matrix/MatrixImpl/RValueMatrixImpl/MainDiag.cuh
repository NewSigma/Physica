/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MainDiag<M>> : public device_obj<RValueVector<MainDiag<M>>> {
        using host_obj = MainDiag<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat) : mat(asStruct(mat)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return std::min(getExpr().getCol(), getExpr().getRow()); }
    };

    template<Matrix M>
    __device__ auto device_obj<MainDiag<M>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        return getExpr().calc(index, index, block);
    }
    template<Matrix M>
    __host__ __device__ auto&& device_obj<MainDiag<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<MainDiag<M>>> : public Traits<MainDiag<M>> {
        static_assert(std::is_reference<M>::value);
    };
}
