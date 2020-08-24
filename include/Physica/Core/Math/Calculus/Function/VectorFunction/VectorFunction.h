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
#ifndef PHYSICA_VECTORFUNCTION_H
#define PHYSICA_VECTORFUNCTION_H

#include <vector>
#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunction.h"

namespace Physica::Core {
    /*!
     * A function can be described using either tree or vector.
     *
     * @class VectorFunction is described using vector, which has a good utilization on cache so it is
     * faster than @class TreeFunction.
     */
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class VectorFunction : public AbstractFunction<type, errorTrack> {
        typedef AbstractFunction<type, errorTrack> Base;
        typedef std::vector<FunctionType> FunctionTypeVector;
        typedef std::vector<Scalar<type, errorTrack>*> ValueVector;
        /*!
         * @class Accessor will read data and turn @class TreeFunctionData into vectors in @class VectorFunction.
         */
        class Accessor {
            FunctionTypeVector& typeVector;
            ValueVector& valueVector;
            const TreeFunction<type, errorTrack>& treeFunc;
        public:
            explicit Accessor(VectorFunction& vectorFunc, const TreeFunction<type, errorTrack>& treeFunc);
        private:
            void access(TreeFunctionData<type, errorTrack>* data);
        };
        /*!
         * \typeVector stores @enum FunctionType in the order that we access a @class TreeFunctionData using DFS.
         * (TreeFunctionData is implemented using binary tree.)
         *
         * e.g. Here is a tree, whose left nodes are always on the top.
         *
         *    sin - x
         *  /
         * *      y
         *  \   /
         *   log
         *      \
         *       cos - x
         *
         * We construct a @class VectorFunctionData from @class TreeFunctionData, so we have a \typeVector:
         * *, sin, value, log, value, cos, value
         */
        FunctionTypeVector typeVector;
        /*!
         * \valueVector stores value we will encounter during the access of a @class TreeFunctionData using DFS.
         *
         * e.g.
         * The example is same to typeVector's example and valueVector contains: x, y and x.
         */
        ValueVector valueVector;
        /*!
         * \valueIte is only used in solve(). \valueIte points to the first value that have not been used.
         */
        typename ValueVector::const_iterator valueIte;
    public:
        explicit VectorFunction(const TreeFunction<type, errorTrack>& treeFunc);
        VectorFunction(const VectorFunction& f);
        VectorFunction(VectorFunction&& f) noexcept;
        ~VectorFunction() = default;
        /* Operations */
        Scalar<type, errorTrack> solve() const { valueIte = valueVector.cbegin(); return solveImpl(typeVector.cbegin()); }
    private:
        Scalar<type, errorTrack> solveImpl(FunctionTypeVector::const_iterator& typeIte) const;
        /* Friends */
        friend class Accessor;
    };
}

#include "VectorFunctionImpl/VectorFunctionImpl.h"

#endif
