/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<Conjugate<M>> : public device_obj<RValueMatrix<Conjugate<M>>> {
        using host_obj = Conjugate<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] __host__ __device__ decltype(auto) transpose(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& conjugate(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) hermite(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getRow() const noexcept { return conjugate().getRow(); }
        [[nodiscard]] __host__ __device__ constexpr size_t getCol() const noexcept { return conjugate().getCol(); }
        [[nodiscard]] __host__ __device__ constexpr size_t getOrder() const noexcept { return conjugate().getOrder(); }
    };

    template<Matrix M>
    __device__ auto device_obj<Conjugate<M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        return conjugate().calc(row, col, block).conjugate();
    }

    template<Matrix M>
    __host__ __device__ decltype(auto) device_obj<Conjugate<M>>::transpose(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).conjugate().hermite();
    }

    template<Matrix M>
    __host__ __device__ auto&& device_obj<Conjugate<M>>::conjugate(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }

    template<Matrix M>
    __host__ __device__ decltype(auto) device_obj<Conjugate<M>>::hermite(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).conjugate().transpose();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<Conjugate<M>>> : public Traits<Conjugate<M>> {
        static_assert(!is_device_obj<M>::value);
    };
}
