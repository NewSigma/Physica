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
    template<class VectorType1, class VectorType2>
    class device_obj<VectorExpr<ExpressionType::Add, VectorType1, VectorType2>>
            : public device_obj<RValueVector<VectorExpr<ExpressionType::Add, VectorType1, VectorType2>>> {
        using host_obj = VectorExpr<ExpressionType::Add, VectorType1, VectorType2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector1 = device_obj<VectorType1>;
        using DeviceVector2 = device_obj<VectorType2>;
    public:
        using typename Base::ScalarType;
    private:
        const DeviceVector1& v1;
        const DeviceVector2& v2;
    public:
        __device__ device_obj(const device_obj<RValueVector<VectorType1>>& v1_, const device_obj<RValueVector<VectorType2>>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1_.getLength() == v2_.getLength());
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return ScalarType(v1.calc(index)) + ScalarType(v2.calc(index));
        }
        [[nodiscard]] __device__ size_t getLength() const { return v1.getLength(); }
    };

    template<class Derived, class OtherDerived>
    [[nodiscard]] __device__ inline auto operator+(
            const device_obj<RValueVector<Derived>>& v1, const device_obj<RValueVector<OtherDerived>>& v2) noexcept {
        return device_obj<VectorExpr<ExpressionType::Add, Derived, OtherDerived>>(v1.getDerived(), v2.getDerived());
    }
}
