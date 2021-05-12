/*
 * Copyright 2021 WeiBo He.
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
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class T, class TVector>
    class ODESolver {
    private:
        TVector x;
        DenseMatrix<T> solution;
        T stepSize;
    public:
        ODESolver(const T& start, const T& end, const T& stepSize_, const TVector& initial);
        /* Operations */
        /**
         * \tparam Function
         * A function object like this
         * TVector func(const T& x, const TVector& y)
         */
        template<class Function>
        void rungeKutta4(Function func);
        template<class Function>
        void verlet(Function func, const TVector& initial1);
        /* Getters */
        const TVector& getX() const noexcept { return x; }
        const DenseMatrix<T>& getSolution() const noexcept { return solution; }
    private:
        template<class Function>
        TVector RungeKuttaDy(size_t step, const TVector& dy_dx, Function func);
    };

    template<class T, class TVector>
    ODESolver<T, TVector>::ODESolver(const T& start, const T& end, const T& stepSize_, const TVector& initial)
            : stepSize(stepSize_) {
        assert(start < end);
        const size_t size = static_cast<size_t>(double((end - start) / stepSize));
        x.reserve(size);
        x.append(start);
        solution.resize(initial.getLength(), size);
        solution[0] = initial;
    }

    template<class T, class TVector>
    template<class Function>
    void ODESolver<T, TVector>::rungeKutta4(Function func) {
        const size_t column_1 = solution.getColumn() - 1;
        for (size_t i = 0; i < column_1; ++i) {
            const T& x_i = x[i];
            TVector dy_dx = func(x_i, solution[i]);
            solution[i + 1] = RungeKuttaDy(i, dy_dx, func);
            x.append(x_i + stepSize);
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.524
     */
    template<class T, class TVector>
    template<class Function>
    TVector ODESolver<T, TVector>::RungeKuttaDy (
            size_t step, const TVector& dy_dx, Function func) {
        const TVector k1 = T(0.5) * stepSize * dy_dx;
        const T temp = x[step] + stepSize * T(0.5);
        const TVector& y = solution[step];
        const TVector k2 = stepSize * func(temp, y + k1);
        const TVector k3 = stepSize * func(temp, y + k2 * T(0.5));
        const TVector k4 = T(0.5) * stepSize * func(x[step] + stepSize, y + k3);
        return y + (k1 + k2 + k3 + k4) / T(3);
    }

    template<class T, class TVector>
    template<class Function>
    void ODESolver<T, TVector>::verlet(Function func, const TVector& initial1) {
        x.append(x[0] + stepSize);
        solution[1] = initial1;
        const T stepSize_2 = square(stepSize);
        const size_t column_1 = solution.getColumn() - 1;
        for (size_t i = 1; i < column_1; ++i) {
            const T& x_i = x[i];
            const TVector& y_i = solution[i];
            solution[i + 1] = -solution[i - 1] + y_i * T(2) + func(x_i, y_i) * stepSize_2;
            x.append(x_i + stepSize);
        }
    }
}