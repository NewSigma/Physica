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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<class T>
    class device_obj<RealMatrix<T>> : public device_obj<RValueMatrix<RealMatrix<T>>> {
        using host_obj = RealMatrix<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& mat;
    public:
    __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).real(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class device_obj<ImagMatrix<T>> : public device_obj<RValueMatrix<ImagMatrix<T>>> {
        using host_obj = ImagMatrix<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& mat;
    public:
    __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).imag(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class device_obj<SquaredNormMatrix<T>> : public device_obj<RValueMatrix<SquaredNormMatrix<T>>> {
        using host_obj = SquaredNormMatrix<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& mat;
    public:
        __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).squaredNorm(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const { return mat.calc(row, col).value().squaredNorm(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class device_obj<NormMatrix<T>> : public device_obj<RValueMatrix<NormMatrix<T>>> {
        using host_obj = NormMatrix<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& mat;
    public:
        __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).norm(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const { return mat.calc(row, col).value().norm(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
    };

    template<class T>
    class device_obj<ValueMatrix<T>> : public device_obj<RValueMatrix<ValueMatrix<T>>> {
        using host_obj = ValueMatrix<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<T>& mat;
    public:
        __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc_value(row, col); }
        [[nodiscard]] __device__ ScalarType calc_value(size_t row, size_t col) const { return calc(row, col); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
    };

    template<class T, int GradOrder>
    class device_obj<GradMatrix<T, GradOrder>> : public device_obj<RValueMatrix<GradMatrix<T, GradOrder>>> {
        using host_obj = GradMatrix<T, GradOrder>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const device_obj<T>& mat;
    public:
        __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).template grad<GradOrder>(); }
        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<RealMatrix<T>>> : public Traits<RealMatrix<T>> {};

    template<class T>
    class Traits<device_obj<ImagMatrix<T>>> : public Traits<ImagMatrix<T>> {};

    template<class T>
    class Traits<device_obj<SquaredNormMatrix<T>>> : public Traits<SquaredNormMatrix<T>> {};

    template<class T>
    class Traits<device_obj<NormMatrix<T>>> : public Traits<NormMatrix<T>> {};

    template<class T>
    class Traits<device_obj<ValueMatrix<T>>> : public Traits<ValueMatrix<T>> {};

    template<class T, int GradOrder>
    class Traits<device_obj<GradMatrix<T, GradOrder>>> : public Traits<GradMatrix<T, GradOrder>> {};
}
