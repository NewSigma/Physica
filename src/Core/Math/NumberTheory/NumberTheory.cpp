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
#include <cassert>
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"
#include "Physica/Core/Exception/NotImplementedException.h"

namespace Physica::Core {
    /**
     * Optimize: Make use of DP and the fact that $C_n^m = C_n^(n - m)$
     */
    Integer bernoulli(const Integer& i) {
        Integer result = 1;
        if(i.isZero())
            return result;
        Integer temp = 0;
        for(; temp < i; ++temp)
            result -= combination(i, temp) * bernoulli(temp) / (i - temp + 1);
        return result;
    }

    Integer GCD::run(const Integer& i1, const Integer& i2, Method method) {
        const Integer* bigger = i1 < i2 ? &i2 : &i1;
        const Integer* smaller = i1 < i2 ? &i1 : &i2;
        Integer x(*bigger), y(*smaller);
        switch (method) {
            /**
             * Implementation of Euclidean algorithm.
             * Different from decreasesTechnique(), it has better performance to big integers.
             */
            case Euclidean:
                while (!y.isZero()) {
                    x = x % y;
                    swap(x, y);
                }
                break;
            /**
             * Implementation of decreases technique.
             * Different from euclideanAlgorithm(),
             * it has better performance to small integers because subtract is faster than division.
             *
             * It is recommended that i1 and i2 are less than INT_MAX(approximately) to avoid overflow.
             */
            case Decreases: {
                int shift = 0; //May overflow here
                while(x.isEven() && y.isEven()) {
                    assert(shift < INT_MAX); //Add a check
                    x >>= 1;
                    y >>= 1;
                    ++shift;
                }

                while(!y.isZero()) {
                    x = x - y;
                    swap(x, y);
                }
                x = x << shift;
                break;
            }
            default:
                throw NotImplementedException();
        }
        return x;
    }
}