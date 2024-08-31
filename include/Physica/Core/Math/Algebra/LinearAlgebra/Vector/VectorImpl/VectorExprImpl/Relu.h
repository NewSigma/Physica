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
    class VectorExpr<ExpressionType::Relu, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Relu, VectorType>> {
        using This = VectorExpr<ExpressionType::Relu, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return relu(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    [[nodiscard]] inline auto relu(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Relu, VectorType>(v);
    }
}
