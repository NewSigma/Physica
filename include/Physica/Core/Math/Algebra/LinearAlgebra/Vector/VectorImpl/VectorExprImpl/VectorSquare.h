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
    class VectorExpr<ExpressionType::Square, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Square, VectorType>> {
        using This = VectorExpr<ExpressionType::Square, VectorType>;
        using Base = RValueVector<This>;
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
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return square(v.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const { return square(v.template packet<AnyPacket>(index)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const { return square(v.template packetPartial<AnyPacket>(index, count)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Square, VectorType> square(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Square, VectorType>(v);
    }
}
