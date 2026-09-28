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

#include "../RValueTensor.cuh"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<class V>
    class device_obj<RealTensor<V>> : public device_obj<RValueTensor<RealTensor<V>>> {
        using host_obj = RealTensor<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return tensor.getDerived().calc(index).real(); }

        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return tensor.getDerived().dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return tensor.getDerived().getShape(); }
    };

    template<class V>
    __host__ __device__ decltype(auto) device_obj<RealTensor<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.tensor.getDerived()).values().reals();
    }

    template<class V>
    class device_obj<ImagTensor<V>> : public device_obj<RValueTensor<ImagTensor<V>>> {
        using host_obj = ImagTensor<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return tensor.getDerived().calc(index).imag(); }

        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return tensor.getDerived().dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return tensor.getDerived().getShape(); }
    };

    template<class V>
    __host__ __device__ decltype(auto) device_obj<ImagTensor<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.tensor.getDerived()).values().imags();
    }

    template<class V>
    class device_obj<SquaredNormTensor<V>> : public device_obj<RValueTensor<SquaredNormTensor<V>>> {
        using host_obj = SquaredNormTensor<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return tensor.getDerived().calc(index).squaredNorm(); }

        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return tensor.getDerived().dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return tensor.getDerived().getShape(); }
    };

    template<class V>
    __host__ __device__ decltype(auto) device_obj<SquaredNormTensor<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.tensor.getDerived()).values().squaredNorms();
    }

    template<class V>
    class device_obj<NormTensor<V>> : public device_obj<RValueTensor<NormTensor<V>>> {
        using host_obj = NormTensor<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return tensor.getDerived().calc(index).norm(); }

        [[nodiscard]] __host__ __device__ decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return tensor.getDerived().dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return tensor.getDerived().getShape(); }
    };

    template<class V>
    __host__ __device__ decltype(auto) device_obj<NormTensor<V>>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.tensor.getDerived()).values().norms();
    }

    template<class V>
    class device_obj<ValueTensor<V>> : public device_obj<RValueTensor<ValueTensor<V>>> {
        using host_obj = ValueTensor<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return tensor.getDerived().calc(index).value(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return tensor.getDerived().dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return tensor.getDerived().getShape(); }
    };

    template<class V, int GradOrder>
    class device_obj<GradTensor<V, GradOrder>> : public device_obj<RValueTensor<GradTensor<V, GradOrder>>> {
        using host_obj = GradTensor<V, GradOrder>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<host_obj>>;
        using Ref = add_device_obj_t<V>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> tensor;
    public:
        __host__ __device__ explicit device_obj(Ref tensor_) : tensor(asStruct(tensor_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return tensor.getDerived().calc(index).template grad<GradOrder>(); }
        [[nodiscard]] __device__ decltype(auto) calc_value(const IndexType& index) const { return calc(index).value(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept { return tensor.getDerived().dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return tensor.getDerived().getShape(); }
    };
}

namespace Physica {
    template<class V>
    class Traits<device_obj<RealTensor<V>>> : public Traits<RealTensor<V>> {};

    template<class V>
    class Traits<device_obj<ImagTensor<V>>> : public Traits<ImagTensor<V>> {};

    template<class V>
    class Traits<device_obj<SquaredNormTensor<V>>> : public Traits<SquaredNormTensor<V>> {};

    template<class V>
    class Traits<device_obj<NormTensor<V>>> : public Traits<NormTensor<V>> {};

    template<class V>
    class Traits<device_obj<ValueTensor<V>>> : public Traits<ValueTensor<V>> {};

    template<class V, int GradOrder>
    class Traits<device_obj<GradTensor<V, GradOrder>>> : public Traits<GradTensor<V, GradOrder>> {};
}
