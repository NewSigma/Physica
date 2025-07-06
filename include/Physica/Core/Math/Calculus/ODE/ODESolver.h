/*
 * Copyright 2021-2025 Weibo He.
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

namespace Physica {
    template<Scalar T, size_t Dim>
    class ODESolver {
    public:
        using VectorType = DenseVector<T, Dim>;
        using SolutionType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Vector, Dim>;
    protected:
        VectorND<T> x;
        SolutionType solution;
        T stepSize;
    public:
        ODESolver(const T& start, const T& end, const T& stepSize_, const VectorType& initial);
        /* Operations */
        void rungeKutta4(std::invocable<T, VectorType> auto fn);
        void verlet(std::invocable<T, T> auto fn, const T& initial1);
        void degenerate_numerov(std::invocable<T> auto fn, const T& tangent);
        /* Getters */
        [[nodiscard]] const VectorND<T>& getX() const noexcept { return x; }
        [[nodiscard]] const SolutionType& getSolution() const noexcept { return solution; }
        [[nodiscard]] size_t getNumStep() const noexcept { return x.getLength(); }
        /* Static members */
        [[nodiscard]] static size_t getNumStep(T start, T end, T stepSize);
        static VectorType rungeKutta4_step(T stepSize, T x, const VectorType& sol, std::invocable<T, VectorType> auto fn);
    };

    template<Scalar T, size_t Dim>
    ODESolver<T, Dim>::ODESolver(const T& start, const T& end, const T& stepSize_, const VectorType& initial)
            : stepSize(stepSize_) {
        assert(start < end);
        const size_t size = getNumStep(start, end, stepSize);
        x.reserve(size);
        x.setLength(size);
        x[0] = start;
        solution.resize(initial.getLength(), size);
        auto col = solution.col(0);
        col = initial;
    }

    template<Scalar T, size_t Dim>
    void ODESolver<T, Dim>::rungeKutta4(std::invocable<T, VectorType> auto fn) {
        using FunctionResult = std::invoke_result<decltype(fn), T, VectorType>::type;
        static_assert(FunctionResult::SizeAtCompile == Dim, "[Possible optimization]: Dimention between ODESolver and functor do not match");
        const size_t col_1 = solution.getCol() - 1;
        for (size_t i = 0; i < col_1; ++i) {
            const T& x_i = x[i];
            solution.asArray()[i + 1] = rungeKutta4_step(stepSize, x_i, solution.asArray()[i], fn);
            x[i + 1] = x_i + stepSize;
        }
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:524
     */
    template<Scalar T, size_t Dim>
    auto ODESolver<T, Dim>::rungeKutta4_step(T stepSize, T x, const VectorType& sol, std::invocable<T, VectorType> auto fn) -> VectorType {
        const VectorType k1 = T(0.5) * stepSize * fn(x, sol);
        const T temp = x + stepSize * T(0.5);
        const VectorType k2 = stepSize * fn(temp, sol + k1);
        const VectorType k3 = stepSize * fn(temp, sol + k2 * T(0.5));
        const VectorType k4 = T(0.5) * stepSize * fn(x + stepSize, sol + k3);
        return sol + (k1 + k2 + k3 + k4) / T(3);
    }
    /**
     * Solve the second order ODE that has form: y''(x) = f(x, y(x)).
     * 
     * A simple and stable method
     * 
     * \param initial1
     * Value of y(h)
     * 
     * Reference:
     * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:572
     */
    template<Scalar T, size_t Dim>
    void ODESolver<T, Dim>::verlet(std::invocable<T, T> auto fn, const T& initial1) {
        x[1] = x[0] + stepSize;
        solution(0, 1) = initial1;
        const T stepSize_2 = square(stepSize);
        const size_t col_1 = solution.getCol() - 1;
        for (size_t i = 1; i < col_1; ++i) {
            const T& x_i = x[i];
            const T& y_i = solution(0, i);
            solution(0, i + 1) = -solution(0, i - 1) + y_i * T(2) + fn(x_i, y_i) * stepSize_2;
            x[i + 1] = x_i + stepSize;
        }
    }
    /**
     * Solve the second order ODE that has form: y''(x) = f(x) * y(x).
     * 
     * Less computational effert and better precision rank than Runge-Kutta4.
     * 
     * \param func
     * The function object of f(x)
     * 
     * \param tangent
     * Tangent value at x[0]
     * 
     * Reference:
     * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:573
     */
    template<Scalar T, size_t Dim>
    void ODESolver<T, Dim>::degenerate_numerov(std::invocable<T> auto fn, const T& tangent) {
        const T x0 = x[0];
        x[1] = x0 + stepSize;
        const T stepSize_2 = square(stepSize);
        const T stepSize_2_12 = stepSize_2 * T(1.0 / 12);
        /* Get y(stepSize) */ {
            const T stepSize_2_6 = stepSize_2 * T(1.0 / 6);
            const T f_minus_step = fn(x0 - stepSize);
            const T temp1 = T(1) - stepSize_2_12 * f_minus_step;
            const T temp2 = T(1) - stepSize_2_6 * f_minus_step;
            const T numerator = (T(2) + T(5.0 / 6) * stepSize_2 * fn(x0)) * temp1 * solution(0, 0) + T(2) * stepSize * tangent * temp2;
            const T f_step = fn(x[1]);
            const T denominator = (T(1) - stepSize_2_12 * f_step) * temp2 + (T(1) - stepSize_2_6 * f_step) * temp1;
            solution(0, 1) = numerator / denominator;
        }

        const size_t col_1 = solution.getCol() - 1;
        T w_i_minus1 = solution(0, 0) * (T(1) - stepSize_2_12 * fn(x[0]));
        T w_i = solution(0, 1) * (T(1) - stepSize_2_12 * fn(x[1]));
        T x_i = x[1];
        for (size_t i = 1; i < col_1; ++i) {
            const T w_i_plus1 = T(2) * w_i + stepSize_2 * fn(x_i) * solution(0, i) - w_i_minus1;
            w_i_minus1 = w_i;
            w_i = w_i_plus1;
            
            x_i += stepSize;
            const T factor = T(1) - stepSize_2_12 * fn(x_i);
            solution(0, i + 1) = w_i / factor;
            x[i + 1] = x_i;
        }
    }

    template<Scalar T, size_t Dim>
    size_t ODESolver<T, Dim>::getNumStep(T start, T end, T stepSize) {
        return static_cast<size_t>(double((end - start) / stepSize));
    }
}