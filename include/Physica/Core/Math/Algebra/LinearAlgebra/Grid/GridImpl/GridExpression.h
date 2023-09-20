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

#include "Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h"

namespace Physica::Core {
    template<ExpressionType type, class T1, class T2 = T1> class GridExpression;

    namespace Internal {
        template<class T> class Traits;

        template<ExpressionType type, class Exp1, class Exp2>
        class Traits<GridExpression<type, Exp1, Exp2>> {
            using ScalarType1 = typename Exp1::ScalarType;
            using RealType = typename ScalarType1::RealType;
            using BinaryScalarType = typename BinaryScalarOpReturnType<ScalarType1, typename Exp2::ScalarType>::Type;
        public:
            using ScalarType = typename std::conditional<type == ExpressionType::Abs, RealType, BinaryScalarType>::type;
        };

        template<ExpressionType type, class Exp, class AnyScalar>
        class Traits<GridExpression<type, Exp, ScalarBase<AnyScalar>>> {
        public:
            using ScalarType = typename BinaryScalarOpReturnType<typename Exp::ScalarType, AnyScalar>::Type;
        };
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class GridType1, class GridType2>
    class GridExpression<ExpressionType::Add, GridType1, GridType2>
            : public RValueGrid<GridExpression<ExpressionType::Add, GridType1, GridType2>> {
        using Base = RValueGrid<GridExpression<ExpressionType::Add, GridType1, GridType2>>;
    public:
        using typename Base::ScalarType;
        using typename Base::Index3D;
    private:
        const GridType1& grid1;
        const GridType2& grid2;
    public:
        GridExpression(const RValueGrid<GridType1>& grid1_, const RValueGrid<GridType2>& grid2_)
                : grid1(grid1_.getDerived()), grid2(grid2_.getDerived()) {
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
    template<class GridType, class AnyScalar>
    class GridExpression<ExpressionType::Mul, GridType, ScalarBase<AnyScalar>>
            : public RValueGrid<GridExpression<ExpressionType::Mul, GridType, ScalarBase<AnyScalar>>> {
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
        using Base = RValueGrid<GridExpression<ExpressionType::Mul, GridType, ScalarBase<AnyScalar>>>;
    public:
        using typename Base::ScalarType;
        using typename Base::Index3D;
    private:
        const GridType& grid;
        const AnyScalar& s;
    public:
        GridExpression(const RValueGrid<GridType>& grid_, const AnyScalar& s_)
                : grid(grid_.getDerived()), s(s_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return ScalarType(grid.calc(index)) * ScalarType(s); }
        [[nodiscard]] size_t getDimX() const noexcept { return grid.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return grid.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return grid.getDimZ(); }
    };

    template<class GridType1, class GridType2>
    class GridExpression<ExpressionType::Mul, GridType1, GridType2>
            : public RValueGrid<GridExpression<ExpressionType::Mul, GridType1, GridType2>> {
        using Base = RValueGrid<GridExpression<ExpressionType::Mul, GridType1, GridType2>>;
    public:
        using typename Base::ScalarType;
        using typename Base::Index3D;
    private:
        const GridType1& grid1;
        const GridType2& grid2;
    public:
        GridExpression(const RValueGrid<GridType1>& grid1_, const RValueGrid<GridType2>& grid2_)
                : grid1(grid1_.getDerived()), grid2(grid2_.getDerived()) {
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
    template<class GridType1, class GridType2>
    inline GridExpression<ExpressionType::Add, GridType1, GridType2> operator+(
            const RValueGrid<GridType1>& g1, const RValueGrid<GridType2>& g2) {
        return GridExpression<ExpressionType::Add, GridType1, GridType2>(g1.getDerived(), g2.getDerived());
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class GridType, class ScalarType>
    inline GridExpression<ExpressionType::Mul, GridType, ScalarBase<ScalarType>> operator*(
            const RValueGrid<GridType>& g, const ScalarBase<ScalarType>& s) {
        return GridExpression<ExpressionType::Mul, GridType, ScalarBase<ScalarType>>(g.getDerived(), s.getDerived());
    }

    template<class ScalarType, class GridType>
    inline GridExpression<ExpressionType::Mul, GridType, ScalarBase<ScalarType>> operator*(
            const ScalarBase<ScalarType>& s, const RValueGrid<GridType>& g) {
        return g * s;
    }

    template<class GridType1, class GridType2>
    inline GridExpression<ExpressionType::Mul, GridType1, GridType2> hadamard(
            const RValueGrid<GridType1>& g1, const RValueGrid<GridType2>& g2) {
        return GridExpression<ExpressionType::Mul, GridType1, GridType2>(g1.getDerived(), g2.getDerived());
    }
}
