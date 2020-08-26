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
#ifndef PHYSICA_VECTORFUNCTIONIMPL_H
#define PHYSICA_VECTORFUNCTIONIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    VectorFunction<type, errorTrack>::Accessor::Accessor(VectorFunction <type, errorTrack> &vectorFunc, const TreeFunction <type, errorTrack> &treeFunc)
            : typeVector(vectorFunc.typeVector), valueVector(vectorFunc.valueVector), treeFunc(treeFunc) {
        access(&treeFunc.getTree());
    }
    /*!
     * access() visits a child tree of @class TreeFunction.
     */
    template<ScalarType type, bool errorTrack>
    void VectorFunction<type, errorTrack>::Accessor::access(TreeFunctionData<type, errorTrack>* data) {
        FunctionType funcType = data->getType();
        typeVector.push_back(funcType);
        if(funcType == Value) {
            const size_t pos = treeFunc.getVariablePos();
            if(pos == 0) {
                pos = treeFunc.getConstantPos();
                valueVector.push_back(&Base::constants[pos]);
            }
            else
                valueVector.push_back(&Base::variables[pos]);
        }
        else {
            access(data->getLeft());
            auto right = data->getRight();
            if(right)
                access(right);
        }
    }

    template<ScalarType type, bool errorTrack>
    VectorFunction<type, errorTrack>::VectorFunction(const TreeFunction<type, errorTrack>& treeFunc) : Base(treeFunc) {
        Accessor(*this, treeFunc);
    }

    template<ScalarType type, bool errorTrack>
    VectorFunction<type, errorTrack>::VectorFunction(const VectorFunction& f)
            : AbstractFunction<type, errorTrack>(f)
            , typeVector(f.typeVector)
            , valueVector(f.valueVector)
            , valueIte(f.valueIte) {}

    template<ScalarType type, bool errorTrack>
    VectorFunction<type, errorTrack>::VectorFunction(VectorFunction&& f) noexcept
            : AbstractFunction<type, errorTrack>(f)
            , typeVector(std::move(f.typeVector))
            , valueVector(std::move(f.valueVector))
            , valueIte(f.valueIte) {}

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> VectorFunction<type, errorTrack>::operator()(const Scalar<type, errorTrack>& s) const {
        AbstractFunction<type, errorTrack>::setVariable(s, 0);
        return solve();
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> VectorFunction<type, errorTrack>::solveImpl(
            typename VectorFunction<type, errorTrack>::FunctionTypeVector::const_iterator& typeIte) const {
        switch(*typeIte++) {
            case Value:
                return *valueIte++;
            case Add:
                return solve(typeIte) + solve(typeIte);
            case Sub: {
                auto s1 = solve(typeIte);
                auto s2 = solve(typeIte);
                return s1 - s2;
            }
            case Mul:
                return solve(typeIte) * solve(typeIte);
            case Div: {
                auto s1 = solve(typeIte);
                auto s2 = solve(typeIte);
                return s1 / s2;
            }
            case Square:
                return square(solve(typeIte));
            case Reciprocal:
                return reciprocal(solve(typeIte));
            case Sqrt:
                return sqrt(solve(typeIte));
            case Factorial:
                return factorial(solve(typeIte));
            case Ln:
                return ln(solve(typeIte));
            case Log:
                return log(solve(typeIte));
            case Exp:
                return exp(solve(typeIte));
            case Cos:
                return cos(solve(typeIte));
            case Sin:
                return sin(solve(typeIte));
            case Tan:
                return tan(solve(typeIte));
            case Sec:
                return sec(solve(typeIte));
            case Csc:
                return csc(solve(typeIte));
            case Cot:
                return cos(solve(typeIte));
            case ArcCos:
                return arccos(solve(typeIte));
            case ArcSin:
                return arcsin(solve(typeIte));
            case ArcTan:
                return arctan(solve(typeIte));
            case ArcSec:
                return arcsec(solve(typeIte));
            case ArcCsc:
                return arccsc(solve(typeIte));
            case ArcCot:
                return arccot(solve(typeIte));
            case Cosh:
                return cosh(solve(typeIte));
            case Sinh:
                return sinh(solve(typeIte));
            case Tanh:
                return tanh(solve(typeIte));
            case Sech:
                return sech(solve(typeIte));
            case Csch:
                return csch(solve(typeIte));
            case Coth:
                return coth(solve(typeIte));
            case ArcCosh:
                return arccosh(solve(typeIte));
            case ArcSinh:
                return arcsinh(solve(typeIte));
            case ArcTanh:
                return arctanh(solve(typeIte));
            case ArcSech:
                return arcsech(solve(typeIte));
            case ArcCsch:
                return arccsch(solve(typeIte));
            case ArcCoth:
                return arccoth(solve(typeIte));
        }
    }
}

#endif
