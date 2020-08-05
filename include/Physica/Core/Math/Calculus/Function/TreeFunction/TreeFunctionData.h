/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TREEFUNCTIONDATA_H
#define PHYSICA_TREEFUNCTIONDATA_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Math/Calculus/Function/FunctionType.h"

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    class TreeFunctionPrinter;

    template<ScalarType type, bool errorTrack>
    class TreeFunction;
    /*!
     * @class TreeFunctionData is the internal data type of a TreeFunction.
     */
    template<ScalarType scalarType = MultiPrecision, bool errorTrack = true>
    class TreeFunctionData {
    private:
        union {
            struct {
                TreeFunctionData* left;
                TreeFunctionData* right;
            };
            struct {
                //value must be allocated by @class Function. value must not be deleted by TreeFunctionData.
                Scalar<scalarType, errorTrack>* value{};
                void* placeHolder{};
            };
        };
        //If type equals to @enum Value, we use the second struct in the union.
        FunctionType type;
    public:
        TreeFunctionData(FunctionType type, TreeFunctionData&& left);
        TreeFunctionData(FunctionType type, TreeFunctionData&& left, TreeFunctionData&& right);
        TreeFunctionData(TreeFunctionData&& func) noexcept;
        ~TreeFunctionData();

        TreeFunctionData& operator=(TreeFunctionData&& f) noexcept;
        /* Getters */
        [[nodiscard]] FunctionType getType() const noexcept { return type; }
        [[nodiscard]] const TreeFunctionData* getLeft() const { return getType() == Value ? nullptr : left; }
        [[nodiscard]] const TreeFunctionData* getRight() const { return getType() == Value ? nullptr : right; }
        [[nodiscard]] const Scalar<scalarType, errorTrack>* getValue() const { return getType() == Value ? value : nullptr; }
    private:
        explicit TreeFunctionData(Scalar<scalarType, errorTrack>* value);
        TreeFunctionData(const TreeFunctionData& func);

        TreeFunctionData& operator=(const TreeFunctionData& func);

        [[nodiscard]] Scalar<scalarType, errorTrack> solve() const;
        friend class TreeFunction<scalarType, errorTrack>;
        friend class TreeFunctionPrinter<scalarType, errorTrack>;
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionDataImpl.h"

#endif