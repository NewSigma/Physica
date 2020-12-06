/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_MESH_H
#define PHYSICA_MESH_H

#include <cstddef>
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/Element.h"
#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica::Core {
    template<int dim>
    class Mesh {
        /**
         * Number of all elements in FEM.
         */
        size_t elementCount;
        /**
         * Number of all nodes in FEM.
         */
        size_t nodeCount;
        /**
         * A array whose length is elementCount, stores elements.
         */
        Element<dim>* elements;
        Point<dim>* nodes;
    public:
        Mesh();
        Mesh(const Mesh<dim>&) = delete;
        Mesh(Mesh<dim>&&) noexcept = delete;
        ~Mesh();
        /* Operators */
        Mesh& operator=(const Mesh<dim>&) = delete;
        Mesh& operator=(Mesh<dim>&&) = delete;
    };
}

#include "MeshImpl.h"

#endif
