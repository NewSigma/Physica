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

#include <cassert>
#include <functional>

namespace Physica::Core {
    template<class ScalarType> class ScalarBase;
    template<class VectorType> class VectorBase;
}

namespace Physica::Core::Math {
    template<class ScalarType, class Function, class VectorType>
    class ConjugateGradient {
        Function func;
        VectorType x;
        ScalarType epsilon;
        ScalarType minStepSize;
        ScalarType minimal;
    public:
        ConjugateGradient(Function func_,
                          const VectorBase<VectorType>& x_,
                          const ScalarBase<ScalarType>& epsilon_,
                          const ScalarBase<ScalarType>& minStepSize_);
        ~ConjugateGradient() = default;
        /* Operations */
        ScalarType compute();
    private:
        [[nodiscard]] VectorType gradient();
        [[nodiscard]] ScalarType climbHill(const VectorType& direction);
    };

    template<class ScalarType, class Function, class VectorType>
    ConjugateGradient<ScalarType, Function, VectorType>::ConjugateGradient(Function func_,
                                                                           const VectorBase<VectorType>& x_,
                                                                           const ScalarBase<ScalarType>& epsilon_,
                                                                           const ScalarBase<ScalarType>& minStepSize_)
            : func(func_), x(x_.getDerived()), epsilon(epsilon_.getDerived()), minStepSize(minStepSize_.getDerived()) {
        assert(minStepSize.isPositive());
    }

    template<class ScalarType, class Function, class VectorType>
    ScalarType ConjugateGradient<ScalarType, Function, VectorType>::compute() {
        const auto length = x.getLength();

        ScalarType y_temp = func(std::cref(x));
        VectorType gradient_temp = gradient();
        VectorType search_direction = -gradient_temp;
        for (size_t i = 0; i <= length; ++i) {
            const ScalarType factor = climbHill(search_direction);
            x += factor * search_direction;

            ScalarType y_new = func(std::cref(x));
            const ScalarType delta = abs(ScalarType(y_temp - y_new));
            y_temp = std::move(y_new);

            if (delta < epsilon)
                break;

            VectorType gradient_new = gradient();
            if (i == length) {
                i = 0;
                search_direction = -gradient_new;
                continue;
            }

            const ScalarType norm = search_direction.norm();
            const ScalarType alpha = gradient_new.squaredNorm() / gradient_temp.squaredNorm();
            search_direction = norm * alpha * search_direction - gradient_new;
            gradient_temp = std::move(gradient_new);

            if (!(gradient_temp * search_direction).isNegative()) {
                i = 0;
                search_direction = -gradient_temp;
                continue;
            }
        }
        return y_temp;
    }

    template<class ScalarType, class Function, class VectorType>
    VectorType ConjugateGradient<ScalarType, Function, VectorType>::gradient() {
        const size_t length = x.getLength();
        VectorType result(length);
        const ScalarType y0 = func(std::cref(x));
        for (size_t i = 0; i < length; ++i) {
            ScalarType new_x = x[i] + minStepSize;
            std::swap(x[i], new_x);
            const ScalarType y1 = func(std::cref(x));
            result[i] = (y1 - y0) / minStepSize;
            std::swap(x[i], new_x);
        }
        return result;
    }

    template<class ScalarType, class Function, class VectorType>
    ScalarType ConjugateGradient<ScalarType, Function, VectorType>::climbHill(const VectorType& direction) {
        ScalarType y0 = func(std::cref(x));
        ScalarType y;
        ScalarType result = minStepSize;
        bool stop;
        do {
            y = func(x + result * direction);
            stop = y > y0;
            result += minStepSize;
            y0 = std::move(y);
        } while(!stop);
        return result;
    }
}
