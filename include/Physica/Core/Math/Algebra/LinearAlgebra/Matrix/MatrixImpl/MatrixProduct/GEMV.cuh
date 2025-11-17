/*
 * Copyright 2024-2025 Weibo He.
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
    template<Matrix T, Vector U>
    class device_obj<GEMV<T, U>> : public device_obj<RValueVector<GEMV<T, U>>> {
        using host_obj = GEMV<T, U>;
        using Base = device_obj<RValueVector<host_obj>>;
        using This = device_obj<host_obj>;
        using DeviceVector = device_obj<U>;
        using DeviceMatrix = device_obj<T>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    private:
        PlainStruct<const std::remove_cvref_t<T>> mat;
        PlainStruct<const std::remove_cvref_t<U>> vec;
    public:
        __host__ __device__ device_obj(T mat_, U vec_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void assign(Vector auto& target) const;

        void reverse(const Vector auto& grad) const noexcept;

        auto values() const noexcept { return mat.getDerived().values() * vec.getDerived().values(); }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return mat.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return vec.getDerived(); }
    };

    template<Matrix T, Vector U>
    __host__ __device__ device_obj<GEMV<T, U>>::device_obj(T mat_, U vec_) : mat(asStruct(mat_)), vec(asStruct(vec_)) {
        assert(getLHS().getCol() == getRHS().getLength());
    }

    template<Matrix T, Vector U>
    __device__ auto device_obj<GEMV<T, U>>::calc(size_t index) const -> ScalarType {
        return getLHS().row(index) * getRHS();
    }

    template<Matrix T, Vector U>
    __host__ __device__ void device_obj<GEMV<T, U>>::assign(Vector auto& target) const {
        if constexpr (IsHost())
            Base::assign(target);
        else {
            if constexpr (MatrixOption::isColMatrix<T>()) {
                Base::assert_assign(target);
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

    template<Matrix T, Vector U>
    void device_obj<GEMV<T, U>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        const auto& m = mat.getDerived();
        const auto& v = vec.getDerived();
        Base::assert_assign(grad);
        if constexpr (ReverseDiff<T>)
            m.reverse(g * v.values().transpose());
        if constexpr (ReverseDiff<U>)
            v.reverse(m.values().transpose() * g);
    }
    // FIXME: Turn it into a member function once we dump to C++23
    template<Matrix T, Vector U>
    [[nodiscard]] __host__ __device__ auto operator*(T&& m, U&& v) noexcept requires(std::remove_cvref_t<T>::RowAtCompile != 1 && CUDA<T> && CUDA<U>) {
        return device_obj<GEMV<T&&, U&&>>(std::forward<T>(m), std::forward<U>(v));
    }
}

namespace Physica {
    template<Matrix T, Vector U>
    class Traits<device_obj<GEMV<T, U>>> : public Traits<GEMV<T, U>> {};
}
