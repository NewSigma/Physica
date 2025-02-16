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

#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearSystem/LinearSystem.h"
#include "NumCharacter.h"

namespace Physica {
    /**
     * Reference:
     * [1] 贾俊平、何晓群、金勇进. 统计学(第六版)[M]. 北京:中国人民大学出版社, 2015:280-282.
     */
    template<Scalar T>
    class LinearFit {
        using VectorType = VectorND<T>;
        using ScalarPair = std::pair<T, T>;
    public:
        [[nodiscard]] static ScalarPair fit(const VectorType& x, const VectorType& y);
        [[nodiscard]] static T relatedCoeff(const VectorType& x, const VectorType& y);
        [[nodiscard]] static ScalarPair deviation(const VectorType& x, const VectorType& y, ScalarPair pair);

        [[nodiscard]] static VectorType polyfit(const VectorType& x, const VectorType& y, int order);
        [[nodiscard]] static VectorType polyval(const VectorType& x, const VectorType& polyCoeff);
    };

    template<Scalar T>
    LinearFit<T>::ScalarPair LinearFit<T>::fit(const VectorType& x, const VectorType& y) {
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
    T LinearFit<T>::relatedCoeff(const VectorType& x, const VectorType& y) {
        return covariance(x, y) / sqrt(x.variance() * y.variance());
    }

    template<Scalar T>
    LinearFit<T>::ScalarPair LinearFit<T>::deviation(
            const VectorType& x, const VectorType& y, ScalarPair pair) {
        const T mean_sse = square(y - pair.first * x - pair.second).sum() / T(x.getLength() - 2);
        const T mean_x = x.mean();
        const T ssx = square(x - mean_x).sum();
        const T deviation_k = sqrt(mean_sse / ssx);
        const T deviation_b = sqrt(mean_sse);
        return std::make_pair(deviation_k, deviation_b);
    }
    /**
     * References:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:496-499
     */
    template<Scalar T>
    LinearFit<T>::VectorType LinearFit<T>::polyfit(const VectorType& x, const VectorType& y, int order) {
        using MatrixType = DenseMatrix<T, MatrixOption::Row | MatrixOption::Vector>;
        assert(x.getLength() == y.getLength() && "[Error]: Dimensions do not match");
        assert(order > 1 && "[Error]: Invalid order");
        const size_t length = x.getLength();
        const int row = order + 1;
        MatrixType working(row, row + 1);
        VectorType x_r(length, 1);
        for (int r = 0; r < row; ++r) {
            working(r, row) = y * x_r;
            VectorType x_c(length, 1);
            for (int c = 0; c < row; ++c) {
                working(r, c) = x_r * x_c;
                x_c = hadamard(x, x_c);
            }
            x_r = hadamard(x, x_r);
        }

        LinearSystem<T> equs(std::move(working));
        equs.gaussJordanPartial();
        return equs.getSolution();
    }

    template<Scalar T>
    LinearFit<T>::VectorType LinearFit<T>::polyval(const VectorType& x, const VectorType& polyCoeff) {
        VectorType result(x.getLength(), 0);
        VectorType x1(x.getLength(), 1);
        for (auto coeff : polyCoeff) {
            result += coeff * x1;
            x1 = hadamard(x1, x);
        }
        return result;
    }
}
