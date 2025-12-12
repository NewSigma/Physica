/*
 * Copyright 2025 Weibo He.
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
    class device_obj<ImagVector<V>> : public RValueVector<ImagVector<V>> {
        using host_obj = ImagVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<V>> v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<V>& v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t s) const { return v.getDerived().calc(s).imag(); }
        [[nodiscard]] __device__ Tv calc_value(size_t s) const { return v.getDerived().calc_value(s).imag(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    class device_obj<SquaredNormVector<V>> : public device_obj<RValueVector<SquaredNormVector<V>>> {
        using host_obj = SquaredNormVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<V>> v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<V>& v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t s) const { return v.getDerived().calc(s).squaredNorm(); }
        [[nodiscard]] __device__ Tv calc_value(size_t s) const { return v.getDerived().calc_value(s).squaredNorm(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    class device_obj<NormVector<V>> : public device_obj<RValueVector<NormVector<V>>> {
        using host_obj = NormVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<V>> v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<V>& v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t s) const { return v.getDerived().calc(s).norm(); }
        [[nodiscard]] __device__ Tv calc_value(size_t s) const { return v.getDerived().calc_value(s).norm(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V>
    class device_obj<ValueVector<V>> : public device_obj<RValueVector<ValueVector<V>>> {
        using host_obj = ValueVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<V>> v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<V>& v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t s) const { return v.getDerived().calc_value(s); }
        [[nodiscard]] __device__ T calc_value(size_t s) const { return calc(s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V, int GradOrder>
    class device_obj<GradVector<V, GradOrder>> : public device_obj<RValueVector<GradVector<V, GradOrder>>> {
        using host_obj = GradVector<V, GradOrder>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<V>> v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<V>& v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t s) const { return v.getDerived().calc(s).template grad<GradOrder>(); }
        [[nodiscard]] __device__ Tv calc_value(size_t s) const { return calc(s).value(); }
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
