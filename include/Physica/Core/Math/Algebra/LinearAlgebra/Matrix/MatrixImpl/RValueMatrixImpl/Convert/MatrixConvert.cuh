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

#include "../../RValueMatrix.cuh"

namespace Physica {
    template<class M>
    class device_obj<RealMatrix<M>> : public device_obj<RValueMatrix<RealMatrix<M>>> {
        using host_obj = RealMatrix<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return mat.getDerived().calc(row, col).real(); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<class M>
    class device_obj<ImagMatrix<M>> : public device_obj<RValueMatrix<ImagMatrix<M>>> {
        using host_obj = ImagMatrix<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return mat.getDerived().calc(row, col).imag(); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<class M>
    class device_obj<SquaredNormMatrix<M>> : public device_obj<RValueMatrix<SquaredNormMatrix<M>>> {
        using host_obj = SquaredNormMatrix<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return mat.getDerived().calc(row, col).squaredNorm(); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return mat.getDerived().calc(row, col).value().squaredNorm(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<class M>
    class device_obj<NormMatrix<M>> : public device_obj<RValueMatrix<NormMatrix<M>>> {
        using host_obj = NormMatrix<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return mat.getDerived().calc(row, col).norm(); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return mat.getDerived().calc(row, col).value().norm(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<class M>
    class device_obj<ValueMatrix<M>> : public device_obj<RValueMatrix<ValueMatrix<M>>> {
        using host_obj = ValueMatrix<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return mat.getDerived().calc(row, col).value(); }
        [[nodiscard]] __device__ T calc_value(size_t row, size_t col) const { return calc(row, col); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<class M, int GradOrder>
    class device_obj<GradMatrix<M, GradOrder>> : public device_obj<RValueMatrix<GradMatrix<M, GradOrder>>> {
        using host_obj = GradMatrix<M, GradOrder>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat_) : mat(asStruct(mat_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return mat.getDerived().calc(row, col).template grad<GradOrder>(); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };
}

namespace Physica {
    template<class M>
    class Traits<device_obj<RealMatrix<M>>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<device_obj<ImagMatrix<M>>> : public Traits<ImagMatrix<M>> {};

    template<class M>
    class Traits<device_obj<SquaredNormMatrix<M>>> : public Traits<SquaredNormMatrix<M>> {};

    template<class M>
    class Traits<device_obj<NormMatrix<M>>> : public Traits<NormMatrix<M>> {};

    template<class M>
    class Traits<device_obj<ValueMatrix<M>>> : public Traits<ValueMatrix<M>> {};

    template<class M, int GradOrder>
    class Traits<device_obj<GradMatrix<M, GradOrder>>> : public Traits<GradMatrix<M, GradOrder>> {};
}
