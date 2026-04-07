/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"

namespace Physica {
    /**
     * Reference:
     * [1] 贾俊平、何晓群、金勇进. 统计学(第六版)[M]. 北京:中国人民大学出版社, 2015:280-282.
     */
    template<Scalar T>
    class LinearFit {
    public:
        [[nodiscard]] static Vector2D<T> fit(const VectorND<T>& x, const VectorND<T>& y);
        [[nodiscard]] static T relatedCoeff(const VectorND<T>& x, const VectorND<T>& y);
        [[nodiscard]] static Vector2D<T> deviation(const VectorND<T>& x, const VectorND<T>& y, Vector2D<T> params);

        [[nodiscard]] static VectorND<T> polyfit(const VectorND<T>& x, const VectorND<T>& y, int order);
        [[nodiscard]] static VectorND<T> polyval(const VectorND<T>& x, const VectorND<T>& polyCoeff);
    };

    template<Scalar T>
    auto LinearFit<T>::fit(const VectorND<T>& x, const VectorND<T>& y) -> Vector2D<T> {
        assert(x.getLength() == y.getLength());
        const size_t length = x.getLength();
        T xy_sum(0);
        T x_square_sum(0);
        T x_sum(0);
        T y_sum(0);
        for (size_t i = 0; i < length; ++i) {
            const auto& x_i = x[i];
            const auto& y_i = y[i];
            xy_sum += x_i * y_i;
            y_sum += y_i;
            x_square_sum += square(x_i);
            x_sum += x_i;
        }
        const T length_1 = reciprocal(T(length));
        const T numerator = xy_sum - x_sum * y_sum * length_1;
        const T denominator = x_square_sum - square(x_sum) * length_1;
        const T k = numerator / denominator;
        return {k, (y_sum - k * x_sum) * length_1};
    }

    template<Scalar T>
    T LinearFit<T>::relatedCoeff(const VectorND<T>& x, const VectorND<T>& y) {
        return covariance(x, y) / sqrt(x.variance() * y.variance());
    }

    template<Scalar T>
    auto LinearFit<T>::deviation(const VectorND<T>& x, const VectorND<T>& y, Vector2D<T> params) -> Vector2D<T> {
        const T mean_sse = square(y - params[0] * x - params[1]).sum() / T(x.getLength() - 2);
        const T mean_x = x.mean();
        const T ssx = square(x - mean_x).sum();
        const T deviation_k = sqrt(mean_sse / ssx);
        const T deviation_b = sqrt(mean_sse);
        return {deviation_k, deviation_b};
    }
    /**
     * References:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:496-499
     */
    template<Scalar T>
    auto LinearFit<T>::polyfit(const VectorND<T>& x, const VectorND<T>& y, int order) -> VectorND<T> {
        assert(x.getLength() == y.getLength() && "[Error]: Dimensions do not match");
        assert(order > 1 && "[Error]: Invalid order");
        const size_t length = x.getLength();
        const int row = order + 1;
        DenseMatrix<T, MatrixMajor::Row> working(row, row + 1);
        VectorND<T> x_r(length, 1);
        for (int r = 0; r < row; ++r) {
            working[r, row] = y * x_r;
            VectorND<T> x_c(length, 1);
            for (int c = 0; c < row; ++c) {
                working[r, c] = x_r * x_c;
                x_c = hadamard(x, x_c);
            }
            x_r = hadamard(x, x_r);
        }
        return DenseLU<T, true>(working.leftCols(row)).inv() * working.col(row);
    }

    template<Scalar T>
    auto LinearFit<T>::polyval(const VectorND<T>& x, const VectorND<T>& polyCoeff) -> VectorND<T> {
        VectorND<T> result(x.getLength(), 0);
        VectorND<T> x1(x.getLength(), 1);
        for (auto coeff : polyCoeff) {
            result += coeff * x1;
            x1 = hadamard(x1, x);
        }
        return result;
    }
}
