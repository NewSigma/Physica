/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/Scalar/ExprType.h"

namespace Physica::Core {
    template<ExprType type, Grid T1, class T2 = T1> class GridExpression;
    //////////////////////////////////////Add//////////////////////////////////////
    template<Grid T1, Grid T2>
    class GridExpression<ExprType::Add, T1, T2>
            : public RValueGrid<GridExpression<ExprType::Add, T1, T2>> {
        using Base = RValueGrid<GridExpression<ExprType::Add, T1, T2>>;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    private:
        const T1& grid1;
        const T2& grid2;
    public:
        GridExpression(const T1& grid1_, const T2& grid2_) : grid1(grid1_), grid2(grid2_) {
            for (int i = 0; i < 3; ++i)
                assert(grid1.getDim()[i] == grid2.getDim()[i]);
        }
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return ScalarType(grid1.calc(index)) + ScalarType(grid2.calc(index)); }
        [[nodiscard]] size_t getDimX() const noexcept { return grid1.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return grid1.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return grid1.getDimZ(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<Grid T, Scalar U>
    class GridExpression<ExprType::Mul, T, U>
            : public RValueGrid<GridExpression<ExprType::Mul, T, U>> {
        using Base = RValueGrid<GridExpression<ExprType::Mul, T, U>>;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    private:
        const T& grid;
        const U& s;
    public:
        GridExpression(const T& grid_, const U& s_)
                : grid(grid_), s(s_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return ScalarType(grid.calc(index)) * ScalarType(s); }
        [[nodiscard]] size_t getDimX() const noexcept { return grid.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return grid.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return grid.getDimZ(); }
    };

    template<Grid T1, Grid T2>
    class GridExpression<ExprType::Mul, T1, T2>
            : public RValueGrid<GridExpression<ExprType::Mul, T1, T2>> {
        using Base = RValueGrid<GridExpression<ExprType::Mul, T1, T2>>;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    private:
        const T1& grid1;
        const T2& grid2;
    public:
        GridExpression(const T1& grid1_, const T2& grid2_)
                : grid1(grid1_), grid2(grid2_) {
            for (int i = 0; i < 3; ++i)
                assert(grid1.getDim()[i] == grid2.getDim()[i]);
        }
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return ScalarType(grid1.calc(index)) * ScalarType(grid2.calc(index)); }
        [[nodiscard]] size_t getDimX() const noexcept { return grid1.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return grid1.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return grid1.getDimZ(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Add//////////////////////////////////////
    template<Grid T1, Grid T2>
    inline GridExpression<ExprType::Add, T1, T2> operator+(const T1& g1, const T2& g2) {
        return GridExpression<ExprType::Add, T1, T2>(g1, g2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<Grid T, Scalar U>
    inline auto operator*(const T& g, const U& x) {
        return GridExpression<ExprType::Mul, T, U>(g, x);
    }

    template<Grid T, Scalar U>
    inline auto operator*(const U& s, const T& g) {
        return g * s;
    }

    template<Grid T1, Grid T2>
    inline auto hadamard(const T1& g1, const T2& g2) {
        return GridExpression<ExprType::Mul, T1, T2>(g1, g2);
    }
}

namespace Physica {
    template<ExprType type, Core::Grid Exp1, Core::Grid Exp2>
    class Traits<GridExpression<type, Exp1, Exp2>> {
        using ScalarType1 = Exp1::ScalarType;
        using RealType = ScalarType1::RealType;
        using BinaryScalarType = Core::Internal::BinaryScalarOpRtnTy<ScalarType1, typename Exp2::ScalarType>::Type;
    public:
        using ScalarType = std::conditional<type == ExprType::Abs, RealType, BinaryScalarType>::type;
    };

    template<ExprType type, Core::Grid Exp, Core::Scalar U>
    class Traits<GridExpression<type, Exp, U>> {
    public:
        using ScalarType = Core::Internal::BinaryScalarOpRtnTy<typename Exp::ScalarType, U>::Type;
    };
}
