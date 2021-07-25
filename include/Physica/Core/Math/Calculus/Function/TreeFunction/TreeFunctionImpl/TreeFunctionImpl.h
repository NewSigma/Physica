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
#ifndef PHYSICA_TREEFUNCTIONIMPL_H
#define PHYSICA_TREEFUNCTIONIMPL_H

namespace Physica::Core {
    /*!
     * Set the variable count and the constant count of this function, but the exact function form has not been initialized.
     * Use setTree() to initialize the exact form before calling solve() or a exception will be thrown.
     *
     * Optimize: Initialize the tree in the constructor to avoid exceptions.
     */
    template<ScalarOption option, bool errorTrack>
    TreeFunction<option, errorTrack>::TreeFunction(size_t variablesLength, size_t constantsLength)
            : Base(variablesLength, constantsLength) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunction<option, errorTrack>::TreeFunction(const TreeFunction& func)
            : Base(func), tree(func.tree) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunction<option, errorTrack>::TreeFunction(TreeFunction&& func) noexcept
            : Base(std::move(func)), tree(std::move(func.tree)) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunction<option, errorTrack>& TreeFunction<option, errorTrack>::operator=(const TreeFunction<option, errorTrack>& func) {
        if(this != &func) {
            Base::operator=(func);
            tree = func.tree;
        }
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunction<option, errorTrack>& TreeFunction<option, errorTrack>::operator=(TreeFunction<option, errorTrack>&& func) noexcept {
        Base::operator=(std::move(tree));
        tree = std::move(func.tree);
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> TreeFunction<option, errorTrack>::operator()(Scalar<option, errorTrack> s1) const {
        AbstractFunction<option, errorTrack>::setVariable(std::move(s1), 0);
        return solve();
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> TreeFunction<option, errorTrack>::operator()(Scalar<option, errorTrack> s1
            , Scalar<option, errorTrack> s2) const {
        AbstractFunction<option, errorTrack>::setVariable(std::move(s1), 0);
        AbstractFunction<option, errorTrack>::setVariable(std::move(s2), 1);
        return solve();
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> TreeFunction<option, errorTrack>::operator()(Scalar<option, errorTrack> s1
            , Scalar<option, errorTrack> s2, Scalar<option, errorTrack> s3) const {
        AbstractFunction<option, errorTrack>::setVariable(std::move(s1), 0);
        AbstractFunction<option, errorTrack>::setVariable(std::move(s2), 1);
        AbstractFunction<option, errorTrack>::setVariable(std::move(s3), 2);
        return solve();
    }

    template<ScalarOption option, bool errorTrack>
    void TreeFunction<option, errorTrack>::printTree(std::ostream& os) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        TreeFunctionPrinter printer(*this, os);
        printer.print();
    }
}

#endif
