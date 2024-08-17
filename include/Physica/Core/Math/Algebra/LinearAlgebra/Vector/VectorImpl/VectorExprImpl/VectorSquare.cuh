/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    template<class VectorType>
    class device_obj<VectorExpr<ExpressionType::Square, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Square, VectorType>>> {
        using host_obj = VectorExpr<ExpressionType::Square, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        Physica::PlainStruct<const DeviceVector> v;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType>>& v_) : v(asStruct(v_.getDerived())) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t s) const { return square(v.getDerived().calc(s)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getDerived().getLength(); }
    };

    template<class VectorType>
    [[nodiscard]] __host__ __device__ inline auto square(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Square, VectorType>>(v);
    }
}
