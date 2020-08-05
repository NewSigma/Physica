/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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
        TreeFunction(std::shared_ptr<TreeFunctionData<type, errorTrack>> f, size_t variablesLength, size_t constantsLength);
        TreeFunction(const TreeFunction& func);
        TreeFunction(TreeFunction&& func) noexcept;
        ~TreeFunction() = default;
        /* Operators */
        TreeFunction& operator=(const TreeFunction& func);
        TreeFunction& operator=(TreeFunction&& func) noexcept;
        [[nodiscard]] Scalar<type, errorTrack> operator()(Scalar<type, errorTrack> s);
        [[nodiscard]] Scalar<type, errorTrack> operator()(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2);
        [[nodiscard]] Scalar<type, errorTrack> operator()(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2, Scalar<type, errorTrack> s3);
        friend std::ostream& operator<<(std::ostream& os, const TreeFunction& f);
        /* Operations */
        [[nodiscard]] Scalar<type, errorTrack> solve() const { return tree->solve(); }
        /* Getters */
        void printTree(std::ostream& os);
        [[nodiscard]] const TreeFunctionData<type, errorTrack>& getTree() const { return *tree; }
        [[nodiscard]] inline TreeFunctionData<type, errorTrack> getVariableNode(size_t index) const;
        [[nodiscard]] inline TreeFunctionData<type, errorTrack> getConstantNode(size_t index) const;
    };

    template<ScalarType type, bool errorTrack>
    inline TreeFunctionData<type, errorTrack> TreeFunction<type, errorTrack>::getVariableNode(size_t index) const {
        Q_ASSERT(index < Base::variables.getLength());
        return TreeFunctionData(&Base::variables[index]);
    }

    template<ScalarType type, bool errorTrack>
    inline TreeFunctionData<type, errorTrack> TreeFunction<type, errorTrack>::getConstantNode(size_t index) const {
        Q_ASSERT(index < Base::constants.getLength());
        return TreeFunctionData(&Base::constants[index]);
    }
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionImpl.h"

#endif
