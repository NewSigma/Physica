/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType>
    class PairModelBase final {
    public:
        constexpr static unsigned int Dim = 3;

        using SearchRangeType = Utils::Array<ssize_t, 3>;
        using LatticeMatrix = typename MDCell::LatticeMatrix;
    public:
        template<class Functor>
        static void forParticleInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func);
        [[nodiscard]] static SearchRangeType estimateRange(const MDCell& cell, ScalarType cutoff);
    };

    template<class ScalarType>
    template<class Functor>
    void PairModelBase<ScalarType>::forParticleInRange(const SearchRangeType& range, const LatticeMatrix& lattice, Functor func) {
        using VectorType = Vector<ScalarType, Dim>;
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);

        VectorType v1, v2, v3;
        for (ssize_t x = -range[0]; x <= range[0]; ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (ssize_t z = -range[2]; z <= range[2]; ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3);
                }
            }
        }
    }

    template<class ScalarType>
    typename PairModelBase<ScalarType>::SearchRangeType
    PairModelBase<ScalarType>::estimateRange(const MDCell& cell, ScalarType cutoff) {
        ssize_t max_x, max_y, max_z;
        /* Estimate range */ {
            const ReciprocalCell reciprocal = cell.reciprocal();
            const auto& lattice = reciprocal.getLattice();
            const ScalarType factor = cutoff * (1 / (2 * M_PI));
            max_x = static_cast<ssize_t>(double(factor * lattice.row(0).norm()) + 1);
            max_y = static_cast<ssize_t>(double(factor * lattice.row(1).norm()) + 1);
            max_z = static_cast<ssize_t>(double(factor * lattice.row(2).norm()) + 1);
        }
        return Utils::Array<ssize_t, 3>{max_x, max_y, max_z};
    }
}
