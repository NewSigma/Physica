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
#ifndef PHYSICA_TREEFUNCTIONDATA_H
#define PHYSICA_TREEFUNCTIONDATA_H

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Calculus/Function/FunctionType.h"
#include "Physica/Core/Math/Calculus/Function/AbstractFunction.h"

namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    class TreeFunctionPrinter;

    template<ScalarOption option, bool errorTrack>
    class TreeFunction;
    /*!
     * @class TreeFunctionData is the internal data type of a TreeFunction.
     *
     * Optimize: Enable copy, remove templates
     */
    template<ScalarOption option = MultiPrecision, bool errorTrack = true>
    class TreeFunctionData {
        static_assert(sizeof(long) == sizeof(void*), "This class is contructed with the assumption");
        typedef AbstractFunction<option, errorTrack> Function;
    private:
        union {
            struct {
                TreeFunctionData* left;
                TreeFunctionData* right;
            };
            struct {
                /*
                 * Positive if index of a variable, negative if index of a constant variable
                 * Never equals to 0, the index starts from 1 or -1
                 */
                long index{};
                /*
                 * This place holder is declared to help clear the data of the first union,
                 * that is, we use 'placeHolder = nullptr' rather than 'right = nullptr' to
                 * indicate that we are operating the second struct instead of the first.
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
        /* Operators */
        TreeFunctionData& operator=(TreeFunctionData&& data) noexcept;
        /* Getters */
        [[nodiscard]] FunctionType getType() const noexcept { return type; }
        [[nodiscard]] const TreeFunctionData* getLeft() const { return getType() == Value ? nullptr : left; }
        [[nodiscard]] const TreeFunctionData* getRight() const { return getType() == Value ? nullptr : right; }
        [[nodiscard]] const Scalar<option, errorTrack>* getValue(const Function& func) const;
    private:
        explicit TreeFunctionData(long index_);
        TreeFunctionData(const TreeFunctionData& data);

        TreeFunctionData& operator=(const TreeFunctionData& data);

        [[nodiscard]] Scalar<option, errorTrack> solve(const Function& func) const;
        friend class TreeFunction<option, errorTrack>;
        friend class TreeFunctionPrinter<option, errorTrack>;
    };

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> operator+(TreeFunctionData<option, errorTrack>&& data1, TreeFunctionData<option, errorTrack>&& data2) {
        return TreeFunctionData<option, errorTrack>(Add, std::move(data1), std::move(data2));
    }

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> operator-(TreeFunctionData<option, errorTrack>&& data1, TreeFunctionData<option, errorTrack>&& data2) {
        return TreeFunctionData<option, errorTrack>(Sub, std::move(data1), std::move(data2));
    }

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> operator*(TreeFunctionData<option, errorTrack>&& data1, TreeFunctionData<option, errorTrack>&& data2) {
        return TreeFunctionData<option, errorTrack>(Mul, std::move(data1), std::move(data2));
    }

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> operator/(TreeFunctionData<option, errorTrack>&& data1, TreeFunctionData<option, errorTrack>&& data2) {
        return TreeFunctionData<option, errorTrack>(Div, std::move(data1), std::move(data2));
    }

    template<ScalarOption option, bool errorTrack>
    inline TreeFunctionData<option, errorTrack> sin(TreeFunctionData<option, errorTrack>&& data) {
        return TreeFunctionData<option, errorTrack>(Sin, std::move(data));
    }
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionDataImpl.h"

#endif