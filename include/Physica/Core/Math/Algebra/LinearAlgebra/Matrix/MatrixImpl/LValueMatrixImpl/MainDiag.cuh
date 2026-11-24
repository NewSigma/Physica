/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../LValueMatrix.cuh"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    class device_obj<MainDiag<M>> : public device_obj<LValueVector<MainDiag<M>>> {
        using host_obj = MainDiag<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat) : mat(asStruct(mat)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return std::min(getExpr().getCol(), getExpr().getRow()); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept { return self.getExpr().data_ptr(index, index); }
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    __host__ __device__ auto&& device_obj<MainDiag<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }
}
