/*
 * Copyright 2020-2021 WeiBo He.
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
#ifndef PHYSICA_TREEFUNCTION_H
#define PHYSICA_TREEFUNCTION_H

#include <utility>
#include <memory>
#include "Physica/Core/Math/Calculus/Function/AbstractFunction.h"
#include "TreeFunctionPrinter.h"

namespace Physica::Core {
    /*!
     * A function can be described using either tree or vector.
     *
     * @class TreeFunction is described using a tree, which has more usages than @class VectorFunction
     * but slower than @class VectorFunction.
     */
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class TreeFunction : public AbstractFunction<type, errorTrack> {
        typedef AbstractFunction<type, errorTrack> Base;

        std::shared_ptr<TreeFunctionData<type, errorTrack>> tree;
    public:
        TreeFunction(size_t variablesLength, size_t constantsLength);
        TreeFunction(const TreeFunction& func);
        TreeFunction(TreeFunction&& func) noexcept;
        ~TreeFunction() = default;
        /* Operators */
        TreeFunction& operator=(const TreeFunction& func);
        TreeFunction& operator=(TreeFunction&& func) noexcept;
        [[nodiscard]] Scalar<type, errorTrack> operator()(Scalar<type, errorTrack> s1) const;
        [[nodiscard]] Scalar<type, errorTrack> operator()(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2) const;
        [[nodiscard]] Scalar<type, errorTrack> operator()(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2, Scalar<type, errorTrack> s3) const;
        /* Operations */
        [[nodiscard]] Scalar<type, errorTrack> solve() const { Q_ASSERT(tree.get()); return tree->solve(*this); };
        /* Getters */
        void printTree(std::ostream& os);
        [[nodiscard]] const TreeFunctionData<type, errorTrack>& getTree() const { return *tree; }
        [[nodiscard]] inline TreeFunctionData<type, errorTrack> getVariableNode(long index) const;
        [[nodiscard]] inline TreeFunctionData<type, errorTrack> getConstantNode(long index) const;
        /* Setters */
        void setTree(std::shared_ptr<TreeFunctionData<type, errorTrack>> p) noexcept { tree = std::move(p); }
    };

    template<ScalarType type, bool errorTrack>
    inline TreeFunctionData<type, errorTrack> TreeFunction<type, errorTrack>::getVariableNode(long index) const {
        assert(0 <= index && static_cast<size_t>(index) < Base::variables.getLength());
        return TreeFunctionData<type, errorTrack>(index + 1);
    }

    template<ScalarType type, bool errorTrack>
    inline TreeFunctionData<type, errorTrack> TreeFunction<type, errorTrack>::getConstantNode(long index) const {
        assert(0 <= index && static_cast<size_t>(index) < Base::constants.getLength());
        return TreeFunctionData<type, errorTrack>(-index - 1);
    }

    template<ScalarType type, bool errorTrack>
    inline std::ostream& operator<<(std::ostream& os, const TreeFunction<type, errorTrack>& f) {
        FunctionPrinter(f, os).print();
        return os;
    }
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionImpl.h"

#endif
