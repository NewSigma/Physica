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
    template<class ScalarType, class PosScalarType, unsigned int Dim> class FreeModel;

    namespace Internal {
        template<class T>
        struct is_free_model {
            constexpr static bool value = false;
        };

        template<class ScalarType, class PosScalarType, unsigned int Dim>
        struct is_free_model<FreeModel<ScalarType, PosScalarType, Dim>> {
            constexpr static bool value = true;
        };
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    class FreeModel final {
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
    public:
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
    };
}
