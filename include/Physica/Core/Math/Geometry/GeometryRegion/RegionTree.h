/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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
