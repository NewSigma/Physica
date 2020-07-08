/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/Function/Function.h"
#include "Physica/Core/Math/Calculus/Function/FunctionPrinter.h"
#include "Physica/Core/Math/Calculus/Function/FunctionTreePrinter.h"
#include "Physica/Core/Math/Calculus/Function/FunctionDerivative.h"

namespace Physica::Core {
    Function::Function(size_t variablesLength, size_t constantsLength)
            : tree(nullptr)
            , variables(new MultiScalar[variablesLength]), constants(new MultiScalar[constantsLength])
            , variablesLength(variablesLength), constantsLength(constantsLength) {}

    Function::Function(FunctionTree&& f, size_t variablesLength, size_t constantsLength)
            : tree(new FunctionTree(std::move(f)))
            , variables(new MultiScalar[variablesLength]), constants(new MultiScalar[constantsLength])
            , variablesLength(variablesLength), constantsLength(constantsLength) {}

    Function::Function(Function&& func) noexcept
            : tree(func.tree), variables(func.variables), constants(func.constants)
            , variablesLength(func.variablesLength), constantsLength(func.constantsLength) {
        func.tree = nullptr;
        func.variables = nullptr;
        func.constants = nullptr;
    }

    Function::~Function() {
        delete tree;
        delete[] variables;
        delete[] constants;
    }

    Function& Function::operator=(Function&& func) noexcept {
        tree = func.tree;
        variables = func.variables;
        func.tree = nullptr;
        func.variables = nullptr;
        variablesLength = func.variablesLength;
        constants = func.constants;
        func.constants = nullptr;
        constantsLength = func.constantsLength;
        return *this;
    }

    MultiScalar Function::operator()(MultiScalar s) {
        setVariable(std::move(s), 0);
        return solve();
    }

    MultiScalar Function::operator()(MultiScalar s1, MultiScalar s2) {
        setVariable(std::move(s1), 0);
        setVariable(std::move(s2), 1);
        return solve();
    }

    MultiScalar Function::operator()(MultiScalar s1, MultiScalar s2, MultiScalar s3) {
        setVariable(std::move(s1), 0);
        setVariable(std::move(s2), 1);
        setVariable(std::move(s3), 2);
        return solve();
    }

    void Function::printTree(std::ostream& os) {
        FunctionTreePrinter(*this, os).print();
    }

    std::ostream& operator<<(std::ostream& os, const Function& f) {
        FunctionPrinter(f, os).print();
        return os;
    }
    /*!
     * Get the variable position from a tree f, starts from 1.
     * return 0 if non-existent.
     */
    size_t Function::getVariablePos(const FunctionTree& f) const {
        for(size_t i = 0; i < variablesLength; ++i) {
            if(variables + i == f.value)
                return i + 1;
        }
        return 0;
    }
}