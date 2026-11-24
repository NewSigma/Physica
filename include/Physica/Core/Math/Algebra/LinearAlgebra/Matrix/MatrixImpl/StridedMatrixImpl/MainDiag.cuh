/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/StridedVector.cuh"
#include "../StridedMatrix.cuh"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    class device_obj<MainDiag<M>> : public device_obj<StridedVector<MainDiag<M>>> {
        using host_obj = MainDiag<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedVector<host_obj>>;
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
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ size_t getStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return std::min(getExpr().getCol(), getExpr().getRow()); }
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ auto&& device_obj<MainDiag<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ auto device_obj<MainDiag<M>>::data_handle(this auto&& self) noexcept {
        return self.getExpr().data_handle();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ size_t device_obj<MainDiag<M>>::getStride() const noexcept {
        return getExpr().getRowStride() + getExpr().getColStride();
    }
}
