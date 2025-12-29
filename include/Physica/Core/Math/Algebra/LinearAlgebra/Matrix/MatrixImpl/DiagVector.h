/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"
#include "RValueMatrix.h"

namespace Physica {
    template<Matrix T, bool isLValueMatrix>
    class DiagVector<T, isLValueMatrix> : public RValueVector<DiagVector<T, isLValueMatrix>> {
        using This = DiagVector<T, isLValueMatrix>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& mat;
    public:
        explicit DiagVector(const T& mat_) : mat(mat_) {}
        DiagVector(const This&) = default;
        DiagVector(This&&) = default;
        ~DiagVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return mat.calc(index, index); }
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
    };

    template<Matrix T>
    class DiagVector<T, true> : public LValueVector<DiagVector<T, true>> {
        using This = DiagVector<T, true>;
        using Base = LValueVector<This>;
    private:
        T& mat;
    public:
        explicit DiagVector(T& mat_) : mat(mat_) {}
        DiagVector(const This&) = default;
        DiagVector(This&&) = default;
        ~DiagVector() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept { return self.mat.data_ptr(index, index); }
    };
}

namespace Physica {
    template<Matrix T, bool isLValueMatrix>
    class Traits<DiagVector<T, isLValueMatrix>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = T::RowAtCompile > T::ColAtCompile ? T::RowAtCompile : T::ColAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
