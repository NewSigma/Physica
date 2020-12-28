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
#ifndef PHYSICA_TREEFUNCTIONDATA_H
#define PHYSICA_TREEFUNCTIONDATA_H

#include "Physica/Core/MultiPrecision/Scalar.h"
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
                //value must be allocated by @class Function and must not be deleted by TreeFunctionData.
                const Scalar<scalarType, errorTrack>* value{};
                /*
                 * This place holder is declared to help clear the data of the first union,
                 * that is, we use 'placeHolder = nullptr' rather than 'right = nullptr' to
                 * indicate thar we are operating the second struct instead of the first.
                 */
                void* placeHolder{};
            };
        };
        //If type equals to @enum Value, we use the second struct in the union.
        FunctionType type;
    public:
        TreeFunctionData(FunctionType type, TreeFunctionData&& left);
        TreeFunctionData(FunctionType type, TreeFunctionData&& left, TreeFunctionData&& right);
        TreeFunctionData(TreeFunctionData&& data) noexcept;
        ~TreeFunctionData();

        TreeFunctionData& operator=(TreeFunctionData&& data) noexcept;
        /* Getters */
        [[nodiscard]] FunctionType getType() const noexcept { return type; }
        [[nodiscard]] const TreeFunctionData* getLeft() const { return getType() == Value ? nullptr : left; }
        [[nodiscard]] const TreeFunctionData* getRight() const { return getType() == Value ? nullptr : right; }
        [[nodiscard]] const Scalar<scalarType, errorTrack>* getValue() const { return getType() == Value ? value : nullptr; }
    private:
        explicit TreeFunctionData(const Scalar<scalarType, errorTrack>* value);
        TreeFunctionData(const TreeFunctionData& data);

        TreeFunctionData& operator=(const TreeFunctionData& data);

        [[nodiscard]] Scalar<scalarType, errorTrack> solve() const;
        friend class TreeFunction<scalarType, errorTrack>;
        friend class TreeFunctionPrinter<scalarType, errorTrack>;
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionDataImpl.h"

#endif