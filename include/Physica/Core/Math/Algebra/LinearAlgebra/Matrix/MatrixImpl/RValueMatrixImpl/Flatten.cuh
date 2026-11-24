/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Flatten.h"
#include "Physica/Core/Parallel/ThreadBlock.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<Flatten<M>> : public device_obj<RValueVector<Flatten<M>>> {
        using host_obj = Flatten<M>;
        using This = device_obj<Flatten<M>>;
        using Base = device_obj<RValueVector<host_obj>>;

        const device_obj<M>& mat;
    protected:
        using typename Base::T;
    public:
        __host__ __device__ device_obj(const device_obj<M>& mat_) : mat(mat_) {}
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
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
    };

    template<Matrix M>
    __device__ auto device_obj<Flatten<M>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor, block);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<Flatten<M>>> : public Traits<Flatten<M>> {};
}
