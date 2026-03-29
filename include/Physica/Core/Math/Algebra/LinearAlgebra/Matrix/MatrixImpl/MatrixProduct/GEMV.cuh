/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/PlainStruct.h"
#include "GEMV.h"

namespace Physica {
    template<Matrix M, Vector V>
    class device_obj<GEMV<M, V>> : public device_obj<RValueVector<GEMV<M, V>>> {
        using host_obj = GEMV<M, V>;
        using Base = device_obj<RValueVector<host_obj>>;
        using This = device_obj<host_obj>;
        using RefM = add_device_obj<M>::type;
        using RefV = add_device_obj<V>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> vec;
    public:
        __host__ __device__ device_obj(RefM mat_, RefV vec_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void assign(Vector auto& target) const;

        [[nodiscard]] __device__ T calc(size_t index) const;
        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M, Vector V>
    __host__ __device__ device_obj<GEMV<M, V>>::device_obj(RefM mat_, RefV vec_) : mat(asStruct(mat_)), vec(asStruct(vec_)) {
        assert(getLHS().getCol() == getRHS().getLength());
    }

    template<Matrix M, Vector V>
    __host__ __device__ void device_obj<GEMV<M, V>>::assign(Vector auto& target) const {
        if constexpr (IsHost())
            Base::assign(target);
        else {
            if constexpr (MatrixMajor::isColMatrix<M>()) {
                target.assert_assign(*this);
                const auto& m = getLHS();
                const auto& v = getRHS();
                target = m.col(0) * v.calc(0);
                for (size_t i = 1; i < v.getLength(); ++i)
                    target += m.col(i) * v.calc(i);
            }
            else
                Base::assign(target);
        }
    }

    template<Matrix M, Vector V>
    __device__ auto device_obj<GEMV<M, V>>::calc(size_t index) const -> T {
        return getLHS().row(index) * getRHS();
    }

    template<Matrix M, Vector V>
    void device_obj<GEMV<M, V>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        const auto& g = grad.values();
        const auto& m = mat.getDerived();
        const auto& v = vec.getDerived();
        Base::assert_assign(grad);
        if constexpr (ReverseDiff<M>)
            m.reverse(g * v.values().transpose());
        if constexpr (ReverseDiff<V>)
            v.reverse(m.values().transpose() * g);
    }

    template<Matrix M, Vector V>
    __host__ __device__ auto device_obj<GEMV<M, V>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V>
    __host__ __device__ auto&& device_obj<GEMV<M, V>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat.getDerived());
    }

    template<Matrix M, Vector V>
    __host__ __device__ auto&& device_obj<GEMV<M, V>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec.getDerived());
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<device_obj<GEMV<M, V>>> : public Traits<GEMV<M, V>> {};
}
