/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTION_H
#define PHYSICA_FUNCTION_H

#include <utility>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"
#include "FunctionTree.h"

namespace Physica::Core {
    class Function {
        FunctionTree* tree;
        MultiScalar* __restrict variables;
        MultiScalar* __restrict constants;
        size_t variablesLength;
        size_t constantsLength;
    public:
        Function(size_t variablesLength, size_t constantsLength);
        Function(FunctionTree&& f, size_t variablesLength, size_t constantsLength);
        Function(const Function& func) = delete;
        Function(Function&& func) noexcept;
        ~Function();
        /* Operators */
        Function& operator=(const Function& func) = delete;
        Function& operator=(Function&& func) noexcept;
        [[nodiscard]] MultiScalar operator()(MultiScalar s);
        [[nodiscard]] MultiScalar operator()(MultiScalar s1, MultiScalar s2);
        [[nodiscard]] MultiScalar operator()(MultiScalar s1, MultiScalar s2, MultiScalar s3);
        friend std::ostream& operator<<(std::ostream& os, const Function& f);
        /* Operations */
        [[nodiscard]] MultiScalar solve() const { return tree->solve(); }
        /* Getters */
        void printTree(std::ostream& os);
        [[nodiscard]] const FunctionTree& getTree() const { return *tree; }
        [[nodiscard]] size_t getVariablesLength() const { return variablesLength; }
        [[nodiscard]] size_t getConstantsLength() const { return constantsLength; }
        [[nodiscard]] inline const MultiScalar& getVariable(size_t index) const;
        [[nodiscard]] inline const MultiScalar& getConstant(size_t index) const;
        [[nodiscard]] inline FunctionTree getVariableNode(size_t index) const;
        [[nodiscard]] inline FunctionTree getConstantNode(size_t index) const;
        [[nodiscard]] size_t getVariablePos (const FunctionTree& f) const;
        /* Setters */
        void setTree(FunctionTree&& f) { if(tree == nullptr) tree = new FunctionTree(std::move(f)); }
        inline void setVariable(MultiScalar s, size_t index);
        inline void setConstant(MultiScalar s, size_t index);
    };

    inline const MultiScalar& Function::getVariable(size_t index) const {
        Q_ASSERT(index < variablesLength);
        return variables[index];
    }

    inline const MultiScalar& Function::getConstant(size_t index) const {
        Q_ASSERT(index < constantsLength);
        return constants[index];
    }

    inline FunctionTree Function::getVariableNode(size_t index) const {
        Q_ASSERT(index < variablesLength);
        return FunctionTree(variables + index);
    }

    inline FunctionTree Function::getConstantNode(size_t index) const {
        Q_ASSERT(index < constantsLength);
        return FunctionTree(constants + index);
    }

    inline void Function::setVariable(MultiScalar s, size_t index) {
        Q_ASSERT(index < variablesLength);
        variables[index] = std::move(s);
    }

    inline void Function::setConstant(MultiScalar s, size_t index) {
        Q_ASSERT(index < constantsLength);
        constants[index] = std::move(s);
    }
}

#endif
