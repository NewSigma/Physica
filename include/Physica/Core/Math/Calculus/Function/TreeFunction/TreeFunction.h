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
    template<ScalarOption option = MultiPrecision, bool errorTrack = true>
    class TreeFunction : public AbstractFunction<option, errorTrack> {
        typedef AbstractFunction<option, errorTrack> Base;

        std::shared_ptr<TreeFunctionData<option, errorTrack>> tree;
    public:
        TreeFunction(size_t variablesLength, size_t constantsLength);
        TreeFunction(const TreeFunction& func);
        TreeFunction(TreeFunction&& func) noexcept;
        ~TreeFunction() = default;
        /* Operators */
        TreeFunction& operator=(const TreeFunction& func);
        TreeFunction& operator=(TreeFunction&& func) noexcept;
        [[nodiscard]] Scalar<option, errorTrack> operator()(Scalar<option, errorTrack> s1) const;
        [[nodiscard]] Scalar<option, errorTrack> operator()(Scalar<option, errorTrack> s1, Scalar<option, errorTrack> s2) const;
        [[nodiscard]] Scalar<option, errorTrack> operator()(Scalar<option, errorTrack> s1, Scalar<option, errorTrack> s2, Scalar<option, errorTrack> s3) const;
        /* Operations */
        [[nodiscard]] Scalar<option, errorTrack> solve() const { Q_ASSERT(tree.get()); return tree->solve(*this); };
        /* Getters */
        void printTree(std::ostream& os);
        [[nodiscard]] const TreeFunction<option, errorTrack>& getTree() const { return *tree; }
        [[nodiscard]] inline TreeFunctionData<option, errorTrack> getVariableNode(long index) const;
        [[nodiscard]] inline TreeFunctionData<option, errorTrack> getConstantNode(long index) const;
        /* Setters */
        void setTree(std::shared_ptr<TreeFunctionData<option, errorTrack>> p) noexcept { tree = std::move(p); }
    };

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> TreeFunction<option, errorTrack>::getVariableNode(long index) const {
        assert(0 <= index && static_cast<size_t>(index) < Base::variables.getLength());
        return TreeFunctionData<option, errorTrack>(index + 1);
    }

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> TreeFunction<option, errorTrack>::getConstantNode(long index) const {
        assert(0 <= index && static_cast<size_t>(index) < Base::constants.getLength());
        return TreeFunctionData<option, errorTrack>(-index - 1);
    }

    template<ScalarOption option, bool errorTrack>
    inline std::ostream& operator<<(std::ostream& os, const TreeFunction<option, errorTrack>& f) {
        FunctionPrinter(f, os).print();
        return os;
    }
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionImpl.h"

#endif
