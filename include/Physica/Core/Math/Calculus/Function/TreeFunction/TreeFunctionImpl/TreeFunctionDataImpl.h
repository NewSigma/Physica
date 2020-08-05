/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TREEFUNCTIONDATAIMPL_H
#define PHYSICA_TREEFUNCTIONDATAIMPL_H

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class TreeFunction.
     * If you want to get a node standing by a value or a variable, ask it from a \class TreeFunction.
     */
    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(Scalar<scalarType, errorTrack>* value) //NOLINT No need to initialize left, right.
            : type(Value), value(value), placeHolder(nullptr) {}

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(FunctionType type, TreeFunctionData&& left)
            : type(type), left(new TreeFunctionData(std::move(left))), right(nullptr) {}

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(FunctionType type, TreeFunctionData&& left, TreeFunctionData&& right)
            : type(type), left(new TreeFunctionData(std::move(left))), right(new TreeFunctionData(std::move(right))) {}
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(const TreeFunctionData &func) : AbstractFunctionData(func.getType()) {
        //Optimize: may be use operator ? and reinterpret_cast to avoid branches.
        if(func.type == Value) {
            value = func.value;
            placeHolder = nullptr;
        }
        else {
            left = new TreeFunctionData(*func.left);
            right = new TreeFunctionData(*func.right);
        }
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::TreeFunctionData(TreeFunctionData&& func) noexcept
            : AbstractFunctionData(func.getType()), left(func.left), right(func.right) {
        func.left = func.right = nullptr;
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>::~TreeFunctionData() {
        if(getType() != Value) {
            delete left;
            delete right;
        }
    }
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>& TreeFunctionData<scalarType, errorTrack>::operator=(const TreeFunctionData& f) {
        if(this != &f) {
            this->~TreeFunctionData();
            if(f.getType() == Value) {
                value = f.value;
                placeHolder = nullptr;
            }
            else {
                left = new TreeFunctionData(*f.left);
                right = new TreeFunctionData(*f.right);
            }
            AbstractFunctionData::operator=(f);
        }
        return *this;
    }

    template<ScalarType scalarType, bool errorTrack>
    TreeFunctionData<scalarType, errorTrack>& TreeFunctionData<scalarType, errorTrack>::operator=(TreeFunctionData&& f) noexcept {
        AbstractFunctionData::operator=(f);
        left = f.left;
        right = f.right;
        f.left = f.right = nullptr;
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
        }
    }
}

#endif
