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

#include "Physica/Utils/Container/Array/Array.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"

namespace Physica::Core {
    class GridBase {
    public:
        using Index3D = Utils::Array<size_t, 3>;
        /* Static members */
        template<class ScalarType, bool IsUnitLattice, class Functor>
        static void forPointInGrid(Index3D dim, const typename PeriodicCell<ScalarType, 3>::LatticeMatrix& lattice, Functor func);
        template<class ScalarType, bool IsUnitLattice, class Functor>
        static void forPointIndexInGrid(Index3D dim, const typename PeriodicCell<ScalarType, 3>::LatticeMatrix& lattice, Functor func);
        template<class Functor> static void forIndexInGrid(Index3D dim, Functor func);
    };

    template<class ScalarType, bool IsUnitLattice, class Functor>
    void GridBase::forPointInGrid(
            Index3D dim, const typename PeriodicCell<ScalarType, 3>::LatticeMatrix& lattice, Functor func) {
        using LatticeMatrix = typename PeriodicCell<ScalarType, 3>::LatticeMatrix;
        using VectorType = Vector<ScalarType, 3>;

        LatticeMatrix sub_lattice{};
        auto a1 = sub_lattice.row(0);
        auto a2 = sub_lattice.row(1);
        auto a3 = sub_lattice.row(2);
        if constexpr (IsUnitLattice)
            sub_lattice = lattice;
        else {
            a1 = lattice.row(0).asVector() * reciprocal(ScalarType(dim[0]));
            a2 = lattice.row(1).asVector() * reciprocal(ScalarType(dim[1]));
            a3 = lattice.row(2).asVector() * reciprocal(ScalarType(dim[2]));
        }

        VectorType v1, v2, v3;
        for (size_t x = 0; x < dim[0]; ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (size_t y = 0; y < dim[1]; ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (size_t z = 0; z < dim[2]; ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3);
                }
            }
        }
    }

    template<class ScalarType, bool IsUnitLattice, class Functor>
    void GridBase::forPointIndexInGrid(
            Index3D dim, const typename PeriodicCell<ScalarType, 3>::LatticeMatrix& lattice, Functor func) {
        static_assert(!ScalarType::isComplex, "[Error]: Position in 3D space can not be complex number");
        using LatticeMatrix = typename PeriodicCell<ScalarType, 3>::LatticeMatrix;
        using VectorType = Vector<ScalarType, 3>;

        LatticeMatrix sub_lattice{};
        auto a1 = sub_lattice.row(0);
        auto a2 = sub_lattice.row(1);
        auto a3 = sub_lattice.row(2);
        if constexpr (IsUnitLattice)
            sub_lattice = lattice;
        else {
            a1 = lattice.row(0).asVector() * reciprocal(ScalarType(dim[0]));
            a2 = lattice.row(1).asVector() * reciprocal(ScalarType(dim[1]));
            a3 = lattice.row(2).asVector() * reciprocal(ScalarType(dim[2]));
        }

        VectorType v1, v2, v3;
        for (size_t x = 0; x < dim[0]; ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (size_t y = 0; y < dim[1]; ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (size_t z = 0; z < dim[2]; ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3, Index3D{x, y, z});
                }
            }
        }
    }

    template<class Functor>
    void GridBase::forIndexInGrid(Index3D dim, Functor func) {
        for (size_t x = 0; x < dim[0]; ++x)
            for (size_t y = 0; y < dim[1]; ++y)
                for (size_t z = 0; z < dim[2]; ++z)
                    func(Index3D{x, y, z});
    }
}
