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
#ifndef PHYSICA_LINEARPROGRAMMING_H
#define PHYSICA_LINEARPROGRAMMING_H

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    /*!
     * Solve the linear programming. Get the maximum of the target function.
     * Pass a minus target function to get the minimum of the function.
     *
     * Reference:
     * [1] Introduction to Algorithms
     */
    class LinearProgramming {
    public:
        enum RestrainType {
            Equal,
            GreaterEqual,
            LessEqual
        };

        enum LPState {
            Normal,
            Infeasiable,
            Infinity
        };
    private:
        /*
         * target[0] is the constant term, target[n] is the coefficient before the n.th variable.
         */
        Vector<> target;
        Matrix<MultiScalar, Row> data;
        //Refactor: devide into two arrays
        size_t* order;
        LPState state;
    public:
        LinearProgramming();
        LinearProgramming(const LinearProgramming& lp) = default;
        LinearProgramming(LinearProgramming&& lp) noexcept;
        ~LinearProgramming();
        /* Operators */
        LinearProgramming& operator=(const LinearProgramming& lp);
        LinearProgramming& operator=(LinearProgramming&& lp) noexcept;
        /* Operations */
        bool addRestrain(Vector<> v, RestrainType type);
        void forceNonNegative(size_t index);
        void solve();
        /* Getters */
        [[nodiscard]] const Vector<>& getTarget() const { return target; }
        [[nodiscard]] const Matrix<MultiScalar, Row>& getData() const { return data; }
        [[nodiscard]] const size_t* getOrder() const { return order; }
        [[nodiscard]] size_t getOrderLength() const { return data.getRow() + data.getColumn(); }
        [[nodiscard]] LPState getState() const { return state; }
        /* Setters */
        bool setTarget(const Vector<>& v);
        bool setTarget(Vector<>&& v);
    private:
    public:
        void initialize();
        void pivot(size_t basic, size_t nonBasic);
        void solveImpl();
        size_t findMinConst() const;
        size_t findPositiveVariable() const;
    };
}

#endif
