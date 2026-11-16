/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.cuh"

namespace Physica {
    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    class device_obj<ImagMatrix<M>> : public device_obj<StridedMatrix<ImagMatrix<M>>> {
        using host_obj = ImagMatrix<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getOrder() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getColStride() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept;
    };

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ device_obj<ImagMatrix<M>>::device_obj(Ref mat_) : mat(asStruct(mat_)) {}

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ void device_obj<ImagMatrix<M>>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ auto device_obj<ImagMatrix<M>>::data_handle(this auto&& self) noexcept {
        return self.mat.getDerived().data_handle()->imag_ptr();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ constexpr size_t device_obj<ImagMatrix<M>>::getRow() const noexcept {
        return mat.getDerived().getRow();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ constexpr size_t device_obj<ImagMatrix<M>>::getCol() const noexcept {
        return mat.getDerived().getCol();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ constexpr size_t device_obj<ImagMatrix<M>>::getOrder() const noexcept {
        return mat.getDerived().getOrder();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ constexpr size_t device_obj<ImagMatrix<M>>::getRowStride() const noexcept {
        return 2 * mat.getDerived().getRowStride();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ constexpr size_t device_obj<ImagMatrix<M>>::getColStride() const noexcept {
        return 2 * mat.getDerived().getColStride();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t device_obj<ImagMatrix<M>>::getRowStrideAtCompile() noexcept {
        return host_obj::getRowStrideAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t device_obj<ImagMatrix<M>>::getColStrideAtCompile() noexcept {
        return host_obj::getColStrideAtCompile();
    }
}
