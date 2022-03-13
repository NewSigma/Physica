/*
 * Copyright 2020-2022 WeiBo He.
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
/*!
 * Bug: if the start of integrate domain is much larger than step size, the result will be 0. May be use taylor series
 * and expend the function to the first order.
 */
namespace Physica::Core {
    //////////////////////////////////Rectangular//////////////////////////////////
    template<class ScalarType>
    Integrate<Rectangular, ScalarType, 1>::Integrate(Base range, ScalarType stepSize_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)) {}

    template<class ScalarType>
    template<class Function>
    ScalarType Integrate<Rectangular, ScalarType, 1>::solve(Function func) const {
        ScalarType result = 0;
        ScalarType start(Base::from()[0]);
        while(start < Base::to()[0]) {
            result += func(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Ladder//////////////////////////////////
    template<class ScalarType>
    Integrate<Ladder, ScalarType, 1>::Integrate(Base range, ScalarType stepSize_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)) {}

    template<class ScalarType>
    template<class Function>
    ScalarType Integrate<Ladder, ScalarType, 1>::solve(Function func) const {
        const ScalarType& from = Base::from()[0];
        const ScalarType& to = Base::to()[0];
        ScalarType result = ((func(from) + func(to)) >> 1);
        ScalarType start(from + stepSize);
        while(start < to) {
            result += func(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Simpson//////////////////////////////////
    template<class ScalarType>
    Integrate<Simpson, ScalarType, 1>::Integrate(Base range, ScalarType stepSize_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)) {}

    template<class ScalarType>
    template<class Function>
    ScalarType Integrate<Simpson, ScalarType, 1>::solve(Function func) const {
        const ScalarType& from = Base::from()[0];
        const ScalarType& to = Base::to()[0];
        ScalarType result = func(from) + func(to);

        ScalarType odd = ScalarType::Zero();
        ScalarType even = ScalarType::Zero();
        bool b = true;
        ScalarType start = from + stepSize;
        while(start < to) {
            ScalarType& toChange = b ? odd : even;
            b = !b;
            toChange += func(start);
            start += stepSize;
        }
        odd <<= 2;
        even <<= 1;
        result += odd + even;
        result *= stepSize;
        result /= ScalarType(3);
        return result;
    }
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    /**
     * \param pointCount Points count on one side of the x-axis. There are same points on both sides.
     */
    template<class ScalarType>
    Integrate<Tanh_Sinh, ScalarType, 1>::Integrate(Base range, ScalarType stepSize_, uint64_t pointCount_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)), pointCount(pointCount_) {}

    template<class ScalarType>
    template<class Function>
    ScalarType Integrate<Tanh_Sinh, ScalarType, 1>::solve(Function func) const {
        const ScalarType& from = Base::from()[0];
        const ScalarType& to = Base::to()[0];

        const ScalarType constant1 = (to - from) >> 1;
        const ScalarType constant2 = constant1 + from;
        const auto& PI_2 = ScalarType(M_PI_2);

        ScalarType result = PI_2 * func(constant2); //Integral value at t = 0
        ScalarType t = 0;
        for(uint64_t i = 0; i < pointCount; ++i) {
            t += stepSize;
            const ScalarType PI_2_sinh = PI_2 * sinh(t);
            const ScalarType cosh_PI_2_sinh = cosh(PI_2_sinh);
            const ScalarType phi = sinh(PI_2_sinh) / cosh_PI_2_sinh;
            const ScalarType phi_derivative = PI_2 * cosh(t) / square(cosh_PI_2_sinh);
            result += phi_derivative * (func(constant2 + constant1 * phi) + func(constant2 - constant1 * phi));
        }
        result *= constant1 * stepSize;
        return result;
    }
    //////////////////////////////////MonteCarlo//////////////////////////////////
    template<class ScalarType, size_t dim>
    Integrate<MonteCarlo, ScalarType, dim>::Integrate(Base range, uint64_t sampleCount_)
            : Base(std::move(range)), sampleCount(sampleCount_) {}

    template<class ScalarType, size_t dim>
    template<class Function, class RandomGenerator>
    ScalarType Integrate<MonteCarlo, ScalarType, dim>::solve(Function func, RandomGenerator& generator) const {
        ScalarType result = 0;
        for (uint64_t i = 0; i < sampleCount; ++i) {
            VectorType x = VectorType::template random<RandomGenerator>(Base::from(), Base::to(), generator);
            result += func(x);
        }

        ScalarType factor = reciprocal(ScalarType(sampleCount));
        for (size_t i = 0; i < Base::from().getLength(); ++i)
            factor *= Base::to()[i] - Base::from()[i];
        result *= factor;
        return result;
    }

    template<class ScalarType, size_t dim>
    template<class Function, class RandomGenerator>
    ScalarType Integrate<MonteCarlo, ScalarType, dim>::solve_e(Function func, RandomGenerator& generator, ScalarType& deviation) const {
        ScalarType total = 0;
        ScalarType total_square = 0;
        for (uint64_t i = 0; i < sampleCount; ++i) {
            const VectorType x = VectorType::template random<RandomGenerator>(Base::from(), Base::to(), generator);
            const ScalarType y = func(x);
            total += y;
            total_square += square(y);
        }

        const ScalarType factor = reciprocal(ScalarType(sampleCount));
        ScalarType variance = total_square * factor - square(total * factor);

        ScalarType factor1 = factor;
        ScalarType factor2 = reciprocal(ScalarType(sampleCount - 1));
        for (size_t i = 0; i < Base::from().getLength(); ++i) {
            const ScalarType delta = Base::to()[i] - Base::from()[i];
            factor1 *= delta;
            factor2 *= square(delta);
        }
        variance *= factor2;
        deviation = sqrt(variance);
        return total * factor1;
    }
}
