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
#ifndef PHYSICA_TREEFUNCTIONDATAIMPL_H
#define PHYSICA_TREEFUNCTIONDATAIMPL_H

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class TreeFunction.
     * If you want to get a node standing by a value or a variable, ask it from a \class TreeFunction.
     */
    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>::TreeFunctionData(long index_) //NOLINT No need to initialize left, right.
            : index(index_), placeHolder(nullptr), type(Value) {}

    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>::TreeFunctionData(FunctionType type, TreeFunctionData&& left)
            : left(new TreeFunctionData(std::move(left))), right(nullptr), type(type) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>::TreeFunctionData(FunctionType type, TreeFunctionData&& left, TreeFunctionData&& right)
            : left(new TreeFunctionData(std::move(left))), right(new TreeFunctionData(std::move(right))), type(type) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
    }
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>::TreeFunctionData(const TreeFunctionData& data) : type(data.type) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        //Optimize: may be use operator ? and reinterpret_cast to avoid branches.
        if(data.type == Value) {
            index = data.index;
            placeHolder = nullptr;
        }
        else {
            left = new TreeFunctionData(*data.left);
            right = new TreeFunctionData(*data.right);
        }
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>::TreeFunctionData(TreeFunctionData&& data) noexcept
            : left(data.left), right(data.right), type(data.type) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        data.left = data.right = nullptr;
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>::~TreeFunctionData() {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        if(getType() != Value) {
            delete left;
            delete right;
        }
    }
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>& TreeFunctionData<option, errorTrack>::operator=(const TreeFunctionData& data) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        if(this != &data) {
            this->~TreeFunctionData();
            type = data.type;
            if(data.type == Value) {
                index = data.index;
                placeHolder = nullptr;
            }
            else {
                left = new TreeFunctionData(*data.left);
                right = new TreeFunctionData(*data.right);
            }
        }
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack>& TreeFunctionData<option, errorTrack>::operator=(TreeFunctionData&& data) noexcept {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        type = data.type;
        left = data.left;
        right = data.right;
        data.left = data.right = nullptr;
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    const Scalar<option, errorTrack>* TreeFunctionData<option, errorTrack>::getValue(const Function& func) const {
        assert(getType() == Value);
        const long normalIndex = std::abs(index) - 1;
        const Array<Scalar<option, errorTrack>>& arr = index > 0 ? func.getVariables() : func.getConstants();
        return &arr[normalIndex];
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> TreeFunctionData<option, errorTrack>::solve(const Function& func) const {
        switch(getType()) {
            case Value:
                return *getValue(func);
            case Add:
                return left->solve(func) + right->solve(func);
            case Sub:
                return left->solve(func) - right->solve(func);
            case Mul:
                return left->solve(func) * right->solve(func);
            case Div:
                return left->solve(func) / right->solve(func);
            case Square:
                return square(left->solve(func));
            case Reciprocal:
                return reciprocal(left->solve(func));
            case Sqrt:
                return sqrt(left->solve(func));
            case Factorial:
                return factorial(left->solve(func));
            case Ln:
                return ln(left->solve(func));
            case Log:
                return log(left->solve(func), right->solve(func));
            case Exp:
                return exp(left->solve(func));
            case Cos:
                return cos(left->solve(func));
            case Sin:
                return sin(left->solve(func));
            case Tan:
                return tan(left->solve(func));
            case Sec:
                return sec(left->solve(func));
            case Csc:
                return csc(left->solve(func));
            case Cot:
                return cot(left->solve(func));
            case ArcCos:
                return arccos(left->solve(func));
            case ArcSin:
                return arcsin(left->solve(func));
            case ArcTan:
                return arctan(left->solve(func));
            case ArcSec:
                return arcsec(left->solve(func));
            case ArcCsc:
                return arccsc(left->solve(func));
            case ArcCot:
                return arccot(left->solve(func));
            case Cosh:
                return cosh(left->solve(func));
            case Sinh:
                return sinh(left->solve(func));
            case Tanh:
                return tanh(left->solve(func));
            case Sech:
                return sech(left->solve(func));
            case Csch:
                return csch(left->solve(func));
            case Coth:
                return coth(left->solve(func));
            case ArcCosh:
                return arccosh(left->solve(func));
            case ArcSinh:
                return arcsinh(left->solve(func));
            case ArcTanh:
                return arctanh(left->solve(func));
            case ArcSech:
                return arcsech(left->solve(func));
            case ArcCsch:
                return arccsch(left->solve(func));
            case ArcCoth:
                return arccosh(left->solve(func));
            default:
                printf("[%s:%d|Fatal]: Not implemented.", __FILENAME__, __LINE__);
                exit(EXIT_FAILURE);
        }
    }
}

#endif
