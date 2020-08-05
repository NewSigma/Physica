/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TREEFUNCTIONIMPL_H
#define PHYSICA_TREEFUNCTIONIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    TreeFunction<type, errorTrack>::TreeFunction(std::shared_ptr<TreeFunctionData<type, errorTrack>> f, size_t variablesLength, size_t constantsLength)
            : Base(variablesLength, constantsLength), tree(std::move(f)) {}

    template<ScalarType type, bool errorTrack>
    TreeFunction<type, errorTrack>::TreeFunction(const TreeFunction& func)
            : Base(func), tree(func.tree) {}

    template<ScalarType type, bool errorTrack>
    TreeFunction<type, errorTrack>::TreeFunction(TreeFunction&& func) noexcept
            : Base(std::move(func)), tree(std::move(func.tree)) {}

    template<ScalarType type, bool errorTrack>
    TreeFunction<type, errorTrack>& TreeFunction<type, errorTrack>::operator=(const TreeFunction<type, errorTrack>& func) {
        if(this != &func) {
            Base::operator=(func);
            tree = func.tree;
        }
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    TreeFunction<type, errorTrack>& TreeFunction<type, errorTrack>::operator=(TreeFunction<type, errorTrack>&& func) noexcept {
        Base::operator=(std::move(tree));
        tree = std::move(func.tree);
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> TreeFunction<type, errorTrack>::operator()(Scalar<type, errorTrack> s) {
        setVariable(std::move(s), 0);
        return solve();
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> TreeFunction<type, errorTrack>::operator()(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2) {
        setVariable(std::move(s1), 0);
        setVariable(std::move(s2), 1);
        return solve();
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> TreeFunction<type, errorTrack>::operator()(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2, Scalar<type, errorTrack> s3) {
        setVariable(std::move(s1), 0);
        setVariable(std::move(s2), 1);
        setVariable(std::move(s3), 2);
        return solve();
    }

    template<ScalarType type, bool errorTrack>
    void TreeFunction<type, errorTrack>::printTree(std::ostream& os) {
        TreeFunctionPrinter printer(*this, os);
        printer.print();
    }

    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const TreeFunction<type, errorTrack>& f) {
        FunctionPrinter(f, os).print();
        return os;
    }
}

#endif
