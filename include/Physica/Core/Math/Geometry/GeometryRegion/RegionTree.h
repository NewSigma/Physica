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
#ifndef PHYSICA_REGIONTREE_H
#define PHYSICA_REGIONTREE_H

#include "GeometryRegion.h"

namespace Physica::Core {
    /*!
     * A tree describes operations between two regions. e.g. and, or ...
     */
    template<int dim>
    class RegionTree : public GeometryRegion<dim> {
        GeometryRegion<dim>* left;
        GeometryRegion<dim>* right;
    public:
        RegionTree(GeometryRegion<dim>&& regionLeft, GeometryRegion<dim>&& regionRight);
        RegionTree(const RegionTree& tree) = delete;
        RegionTree(RegionTree&& tree) noexcept;
        ~RegionTree() override;
        /* Operators */
        RegionTree& operator=(const RegionTree& tree) = delete;
        RegionTree& operator=(RegionTree&& tree) noexcept;
    protected:
        GeometryRegion<dim>* release() override;
    };

    template<int dim>
    RegionTree<dim>::RegionTree(GeometryRegion<dim>&& regionLeft, GeometryRegion<dim>&& regionRight)
            : left(regionLeft.release()), right(regionRight.release()) {}

    template<int dim>
    RegionTree<dim>::RegionTree(RegionTree<dim>&& tree) noexcept : left(tree.left), right(tree.right) {
        tree.left = tree.right = nullptr;
    }

    template<int dim>
    RegionTree<dim>::~RegionTree() {
        delete left;
        delete right;
    }

    template<int dim>
    RegionTree<dim>& RegionTree<dim>::operator=(RegionTree<dim>&& tree) noexcept {
        left = tree.left;
        right = tree.right;
        tree.left = tree.right = nullptr;
        return *this;
    }

    template<int dim>
    GeometryRegion<dim>* RegionTree<dim>::release() {
        return new RegionTree(std::move(*this));
    }
}

#endif
