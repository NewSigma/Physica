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
    template<class T>
    class device_obj<RealVector<T>> : public device_obj<RValueVector<RealVector<T>>> {
        using host_obj = RealVector<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& v_) : v(v_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t s) const { return v.calc(s).real(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t s) const { return v.calc_value(s).real(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class T>
    class device_obj<ImagVector<T>> : public RValueVector<ImagVector<T>> {
        using host_obj = ImagVector<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& v_) : v(v_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t s) const { return v.calc(s).imag(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t s) const { return v.calc_value(s).imag(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class T>
    class device_obj<SquaredNormVector<T>> : public device_obj<RValueVector<SquaredNormVector<T>>> {
        using host_obj = SquaredNormVector<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    private:
        const device_obj<T>& v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& v_) : v(v_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ CoDiff<ScalarType> calc(size_t s) const { return v.calc(s).squaredNorm(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t s) const { return v.calc_value(s).squaredNorm(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class T>
    class device_obj<NormVector<T>> : public device_obj<RValueVector<NormVector<T>>> {
        using host_obj = NormVector<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& v_) : v(v_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ CoDiff<ScalarType> calc(size_t s) const { return v.calc(s).norm(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t s) const { return v.calc_value(s).norm(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class T>
    class device_obj<ValueVector<T>> : public device_obj<RValueVector<ValueVector<T>>> {
        using host_obj = ValueVector<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& v_) : v(v_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t s) const { return v.calc_value(s); }
        [[nodiscard]] __device__ ScalarType calc_value(size_t s) const { return calc(s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class T, int GradOrder>
    class device_obj<GradVector<T, GradOrder>> : public device_obj<RValueVector<GradVector<T, GradOrder>>> {
        using host_obj = GradVector<T, GradOrder>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& v_) : v(v_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t s) const { return v.calc(s).template grad<GradOrder>(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t s) const { return calc(s).value(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<RealVector<T>>> : public Traits<RealVector<T>> {};

    template<class T>
    class Traits<device_obj<ImagVector<T>>> : public Traits<ImagVector<T>> {};

    template<class T>
    class Traits<device_obj<SquaredNormVector<T>>> : public Traits<SquaredNormVector<T>> {};

    template<class T>
    class Traits<device_obj<NormVector<T>>> : public Traits<NormVector<T>> {};

    template<class T>
    class Traits<device_obj<ValueVector<T>>> : public Traits<ValueVector<T>> {};

    template<class T, int GradOrder>
    class Traits<device_obj<GradVector<T, GradOrder>>> : public Traits<GradVector<T, GradOrder>> {};
}
