/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_INTEGRATEIMPL_H
#define PHYSICA_INTEGRATEIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    AbstractIntegrate<type, errorTrack>::AbstractIntegrate(std::shared_ptr<TreeFunction<type, errorTrack>>&& p)
            : p_func(std::move(p)) {}
    ////////////////////////////////////Dynamic////////////////////////////////////
    ////////////////////////////////////1D////////////////////////////////////
    template<ScalarType type, bool errorTrack>
    Integrate<1, type, errorTrack>::Integrate(std::shared_ptr<TreeFunction<type, errorTrack>> p
            , Scalar<type, errorTrack> from, Scalar<type, errorTrack> to)
                : AbstractIntegrate<type, errorTrack>(std::move(p)), from(std::move(from)), to(std::move(to)) {}
    ////////////////////////////////////2D////////////////////////////////////
    ////////////////////////////////////3D////////////////////////////////////
}

#endif