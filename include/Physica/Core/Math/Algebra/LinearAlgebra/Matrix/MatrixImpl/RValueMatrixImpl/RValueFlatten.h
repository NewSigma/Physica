/*
 * Copyright 2022-2024 Weibo He.
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
    template<class MatrixType>
    class RValueFlatten : public RValueVector<RValueFlatten<MatrixType>> {
        using This = RValueFlatten<MatrixType>;

        const MatrixType& mat;
    public:
        using Base = RValueVector<RValueFlatten<MatrixType>>;
        using typename Base::ScalarType;
    public:
        RValueFlatten(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        RValueFlatten(const This&) = delete;
        RValueFlatten(This&&) noexcept = delete;
        ~RValueFlatten() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getColumn(); }
    };

    template<class MatrixType>
    typename RValueFlatten<MatrixType>::ScalarType RValueFlatten<MatrixType>::calc(size_t index) const {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor);
    }
}

namespace Physica {
    template<class MatrixType>
    class Traits<Core::RValueFlatten<MatrixType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile * MatrixType::ColumnAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
