/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.cuh"

namespace Physica {
    template<class V>
    class device_obj<RealVector<V>> : public device_obj<RValueVector<RealVector<V>>> {
        using host_obj = RealVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj_t<V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> v;
    public:
        __host__ __device__ explicit device_obj(Ref v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const { return v.getDerived().calc(i, block).real(); }

        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    __host__ __device__ decltype(auto) device_obj<RealVector<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.v.getDerived()).values().reals();
    }
}

namespace Physica {
    template<class V>
    class Traits<device_obj<RealVector<V>>> : public Traits<RealVector<V>> {};
}
