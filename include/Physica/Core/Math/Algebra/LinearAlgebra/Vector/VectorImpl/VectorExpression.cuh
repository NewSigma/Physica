/*
 * Copyright 2023 WeiBo He.
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
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    class device_obj<VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>>
            : public device_obj<RValueVector<VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>>> {
    public:
        using host_obj = VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>;
    private:
        using Base = device_obj<RValueVector<host_obj>>;
        using DeviceVector = device_obj<VectorType>;
        
        Physica::PlainStruct<DeviceVector> exp;
    public:
        device_obj(const device_obj<RValueVector<VectorType>>& exp_) : exp(asStruct(exp_.getDerived())) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = delete;
        device_obj& operator=(device_obj&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ typename Base::ScalarType calc(size_t index) const {
            return reciprocal(exp.getDerived().calc(index));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const {
            return exp.getDerived().getLength();
        }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    inline device_obj<VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>>
    reciprocal(const device_obj<RValueVector<VectorType>>& v) {
        return v;
    }
}
