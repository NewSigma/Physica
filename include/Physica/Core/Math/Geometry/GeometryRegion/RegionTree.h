/*
 * Copyright 2020-2024 Weibo He.
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
#pragma once

#include "GeometryRegion.h"

namespace Physica::Core {
    /*!
     * A tree describes operations between two regions. e.g. and, or ...
     */
    template<int Dim>
    class RegionTree : public GeometryRegion<Dim> {
        GeometryRegion<Dim>* left;
        GeometryRegion<Dim>* right;
    public:
        RegionTree(GeometryRegion<Dim>&& regionLeft, GeometryRegion<Dim>&& regionRight);
        RegionTree(const RegionTree& tree) = delete;
        RegionTree(RegionTree&& tree) noexcept;
        ~RegionTree() override;
        /* Operators */
        RegionTree& operator=(const RegionTree& tree) = delete;
        RegionTree& operator=(RegionTree&& tree) noexcept;
    protected:
        GeometryRegion<Dim>* release() override;
    };

    template<int Dim>
    RegionTree<Dim>::RegionTree(GeometryRegion<Dim>&& regionLeft, GeometryRegion<Dim>&& regionRight)
            : left(regionLeft.release()), right(regionRight.release()) {}

    template<int Dim>
    RegionTree<Dim>::RegionTree(RegionTree<Dim>&& tree) noexcept : left(tree.left), right(tree.right) {
        tree.left = tree.right = nullptr;
    }

    template<int Dim>
    RegionTree<Dim>::~RegionTree() {
        delete left;
        delete right;
    }

    template<int Dim>
    RegionTree<Dim>& RegionTree<Dim>::operator=(RegionTree<Dim>&& tree) noexcept {
        left = tree.left;
        right = tree.right;
        tree.left = tree.right = nullptr;
        return *this;
    }

    template<int Dim>
    GeometryRegion<Dim>* RegionTree<Dim>::release() {
        return new RegionTree(std::move(*this));
    }
}
