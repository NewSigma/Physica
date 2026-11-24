/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "GeometryRegion.h"

namespace Physica {
    /**
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
