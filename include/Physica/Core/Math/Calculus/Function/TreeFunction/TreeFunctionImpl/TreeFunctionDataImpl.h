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
#ifndef PHYSICA_TREEFUNCTIONDATAIMPL_H
#define PHYSICA_TREEFUNCTIONDATAIMPL_H

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class TreeFunction.
     * If you want to get a node standing by a value or a variable, ask it from a \class TreeFunction.
     */
    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(const Scalar<scalarType, errorTrack>* value) //NOLINT No need to initialize left, right.
            : type(Value), value(value), placeHolder(nullptr) {}

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(FunctionType type, TreeFunctionData&& left)
            : type(type), left(new TreeFunctionData(std::move(left))), right(nullptr) {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(FunctionType type, TreeFunctionData&& left, TreeFunctionData&& right)
            : type(type), left(new TreeFunctionData(std::move(left))), right(new TreeFunctionData(std::move(right))) {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
    }
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(const TreeFunctionData& data) : type(data.type) {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
        //Optimize: may be use operator ? and reinterpret_cast to avoid branches.
        if(data.type == Value) {
            value = data.value;
            placeHolder = nullptr;
        }
        else {
            left = new TreeFunctionData(*data.left);
            right = new TreeFunctionData(*data.right);
        }
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(TreeFunctionData&& data) noexcept
            : type(data.type), left(data.left), right(data.right) {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
        data.left = data.right = nullptr;
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::~TreeFunctionData() {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
        if(getType() != Value) {
            delete left;
            delete right;
        }
    }
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>& TreeFunctionData<scalarType, errorTrack>::operator=(const TreeFunctionData& data) {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
        if(this != &data) {
            this->~TreeFunctionData();
            type = data.type;
            if(data.type == Value) {
                value = data.value;
                placeHolder = nullptr;
            }
            else {
                left = new TreeFunctionData(*data.left);
                right = new TreeFunctionData(*data.right);
            }
        }
        return *this;
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>& TreeFunctionData<scalarType, errorTrack>::operator=(TreeFunctionData&& data) noexcept {
        Q_UNUSED(scalarType)
        Q_UNUSED(errorTrack)
        type = data.type;
        left = data.left;
        right = data.right;
        data.left = data.right = nullptr;
        return *this;
    }

    template<ScalarType scalarType, bool errorTrack>
    Scalar<scalarType, errorTrack> TreeFunctionData<scalarType, errorTrack>::solve() const {
        switch(getType()) {
            case Value:
                return Scalar<scalarType, errorTrack>(*value);
            case Add:
                return left->solve() + right->solve();
            case Sub:
                return left->solve() - right->solve();
            case Mul:
                return left->solve() * right->solve();
            case Div:
                return left->solve() / right->solve();
            case Square:
                return square(left->solve());
            case Reciprocal:
                return reciprocal(left->solve());
            case Sqrt:
                return sqrt(left->solve());
            case Factorial:
                return factorial(left->solve());
            case Ln:
                return ln(left->solve());
            case Log:
                return log(left->solve(), right->solve());
            case Exp:
                return exp(left->solve());
            case Cos:
                return cos(left->solve());
            case Sin:
                return sin(left->solve());
            case Tan:
                return tan(left->solve());
            case Sec:
                return sec(left->solve());
            case Csc:
                return csc(left->solve());
            case Cot:
                return cot(left->solve());
            case ArcCos:
                return arccos(left->solve());
            case ArcSin:
                return arcsin(left->solve());
            case ArcTan:
                return arctan(left->solve());
            case ArcSec:
                return arcsec(left->solve());
            case ArcCsc:
                return arccsc(left->solve());
            case ArcCot:
                return arccot(left->solve());
            case Cosh:
                return cosh(left->solve());
            case Sinh:
                return sinh(left->solve());
            case Tanh:
                return tanh(left->solve());
            case Sech:
                return sech(left->solve());
            case Csch:
                return csch(left->solve());
            case Coth:
                return coth(left->solve());
            case ArcCosh:
                return arccosh(left->solve());
            case ArcSinh:
                return arcsinh(left->solve());
            case ArcTanh:
                return arctanh(left->solve());
            case ArcSech:
                return arcsech(left->solve());
            case ArcCsch:
                return arccsch(left->solve());
            case ArcCoth:
                return arccosh(left->solve());
            default:
                printf("[%s:%d|Fatal]: Not implemented.", __FILENAME__, __LINE__);
                exit(EXIT_FAILURE);
        }
    }
}

#endif
