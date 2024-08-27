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
    class device_obj<VectorExpr<ExpressionType::Minus, VectorType>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Minus, VectorType>>> {
        using host_obj = VectorExpr<ExpressionType::Minus, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
    public:
        using typename Base::ScalarType;
    private:
        const DeviceVector& exp;
    public:
        __device__ device_obj(const device_obj<RValueVector<VectorType>>& exp_) : exp(exp_.getDerived()) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t s) const { return -exp.calc(s); }
        [[nodiscard]] __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class Derived>
    [[nodiscard]] __device__ inline auto operator-(const device_obj<RValueVector<Derived>>& v) noexcept {
        return device_obj<VectorExpr<ExpressionType::Minus, Derived>>(v.getDerived());
    }
}
