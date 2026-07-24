/*
 * Copyright 2025-2026 Weibo He.
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

#include "VectorConvert.h"

namespace Physica {
    template<class V>
    class device_obj<ImagVector<V>> : public device_obj<RValueVector<ImagVector<V>>> {
        using host_obj = ImagVector<V>;
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
        This& operator=(This&&) = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const { return v.getDerived().calc(i, block).imag(); }

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    __host__ __device__ auto device_obj<ImagVector<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.v.getDerived()).values().imags();
    }

    template<class V>
    class device_obj<SquaredNormVector<V>> : public device_obj<RValueVector<SquaredNormVector<V>>> {
        using host_obj = SquaredNormVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using Base::isReverseDiff;
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
        This& operator=(This&&) = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const { return v.getDerived().calc(i, block).squaredNorm(); }

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    __host__ __device__ auto device_obj<SquaredNormVector<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.v.getDerived()).values().squaredNorms();
    }

    template<class V>
    class device_obj<NormVector<V>> : public device_obj<RValueVector<NormVector<V>>> {
        using host_obj = NormVector<V>;
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
        This& operator=(This&&) = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const { return v.getDerived().calc(i, block).norm(); }

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    __host__ __device__ auto device_obj<NormVector<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.v.getDerived()).values().norms();
    }

    template<class V>
    class device_obj<ValueVector<V>> : public device_obj<RValueVector<ValueVector<V>>> {
        using host_obj = ValueVector<V>;
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
        This& operator=(This&&) = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const { return v.getDerived().calc(i, block).value(); }
        [[nodiscard]] __device__ T calc_value(size_t i) const { return calc(i); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V, int GradOrder>
    class device_obj<GradVector<V, GradOrder>> : public device_obj<RValueVector<GradVector<V, GradOrder>>> {
        using host_obj = GradVector<V, GradOrder>;
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
        This& operator=(This&&) = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const { return v.getDerived().calc(i, block).template grad<GradOrder>(); }
        [[nodiscard]] __device__ Tv calc_value(size_t i) const { return calc(i).value(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };
}

namespace Physica {
    template<class V>
    class Traits<device_obj<ImagVector<V>>> : public Traits<ImagVector<V>> {};

    template<class V>
    class Traits<device_obj<SquaredNormVector<V>>> : public Traits<SquaredNormVector<V>> {};

    template<class V>
    class Traits<device_obj<NormVector<V>>> : public Traits<NormVector<V>> {};

    template<class V>
    class Traits<device_obj<ValueVector<V>>> : public Traits<ValueVector<V>> {};

    template<class V, int GradOrder>
    class Traits<device_obj<GradVector<V, GradOrder>>> : public Traits<GradVector<V, GradOrder>> {};
}
