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
#ifndef PHYSICA_GEOMETRYREGION_H
#define PHYSICA_GEOMETRYREGION_H

#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica::Core {
    /*!
     * Describe a region in 2D or 3D. Used tree to implement it.
     */
    template<int dim>
    class GeometryRegion {
        static_assert(dim == 2 || dim == 3, "Dimension of a region must be either 2 or 3.");
    };

    template<>
    class GeometryRegion<2> {
    public:
        enum RegionType {
            And,
            Or,
            Not,
            Xor,
            Square
        };
    private:
        const RegionType type;
    public:
        GeometryRegion(const GeometryRegion& g) = delete;
        GeometryRegion(GeometryRegion&& g) noexcept = delete;
        virtual ~GeometryRegion() = default;
        /* Operators */
        GeometryRegion& operator=(const GeometryRegion& g) = delete;
        GeometryRegion& operator=(GeometryRegion&& g) noexcept = delete;
        /* Getters */
        [[nodiscard]] RegionType getType() const { return type; }
        /*!
         * Return whether this region is a basic region. e.g. a circle, a rectangular etc.
         */
        [[nodiscard]] bool isBasic() const { return type < Xor; }
    protected:
        explicit GeometryRegion(RegionType type);
        /*!
         * New a object and take the ownership of the data in the current object.
         *
         * This function should only be called by \class RegionTree, which assembles regions into its tree.
         */
        virtual GeometryRegion* release() = 0;
    };

    template<>
    class GeometryRegion<3> {
    public:
        enum RegionType {
            And,
            Or,
            Not,
            Xor
        };
    private:
        RegionType type;
    public:
        GeometryRegion(const GeometryRegion& g) = delete;
        GeometryRegion(GeometryRegion&& g) noexcept = delete;
        virtual ~GeometryRegion() = default;
        /* Operators */
        GeometryRegion& operator=(const GeometryRegion& g) = delete;
        GeometryRegion& operator=(GeometryRegion&& g) noexcept = delete;
        /* Getters */
        [[nodiscard]] RegionType getType() const { return type; }
        /*!
         * Return whether this region is a basic region. e.g. a circle, a rectangular etc.
         */
        [[nodiscard]] bool isBasic() const { return type < Xor; }
    protected:
        explicit GeometryRegion(RegionType type);
        /*!
         * New a object and take the ownership of the data in the current object.
         *
         * This function should only be called by \class RegionTree, which assembles regions into its tree.
         */
        virtual GeometryRegion* release() = 0;
    };
}

#endif
