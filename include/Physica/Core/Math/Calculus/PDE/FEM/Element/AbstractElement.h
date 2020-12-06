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
#ifndef PHYSICA_ABSTRACTELEMENT_H
#define PHYSICA_ABSTRACTELEMENT_H

#include <cstddef>
#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica::Core {
    template<size_t dim>
    class AbstractElement {
    public:
        static constexpr size_t dimension = 2;
    private:
        size_t nodesCount;
        /*!
         * A array that stores pointers to nodes, whose length can be accessed from ElementType.
         */
        Point<dim>** nodes;
    public:
        AbstractElement() = delete;
        ~AbstractElement();
        AbstractElement(const AbstractElement&) = delete;
        AbstractElement(AbstractElement&&) noexcept = delete;
        /* Operators */
        AbstractElement& operator=(const AbstractElement&) = delete;
        AbstractElement& operator=(AbstractElement&&) noexcept = delete;
    protected:
        explicit AbstractElement(size_t nodesCount);
    };

    template<size_t dim>
    AbstractElement<dim>::AbstractElement(size_t nodesCount) : nodesCount(nodesCount), nodes(new Point<dim>*[nodesCount]) {}

    template<size_t dim>
    AbstractElement<dim>::~AbstractElement() {
        delete[] nodes;
    }
}

#endif
