/*
 * Copyright 2020-2024 Weibo He.
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

#include "Physica/Core/Math/Statistics/NumCharacter.h"
/*!
 * Bug: if the start of integrate domain is much larger than step size, the result will be 0. May be use taylor series
 * and expend the function to the first order.
 */
namespace Physica::Core {
    //////////////////////////////////Rectangular//////////////////////////////////
    template<Scalar T>
    Integrate<IntegrateMethod::Rectangular, T, 1>::Integrate(Base range, T stepSize_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)) {}

    template<Scalar T>
    template<class Function>
    T Integrate<IntegrateMethod::Rectangular, T, 1>::solve(Function func) const {
        T result = 0;
        T start(Base::from()[0]);
        while(start < Base::to()[0]) {
            result += func(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Ladder//////////////////////////////////
    template<Scalar T>
    Integrate<IntegrateMethod::Ladder, T, 1>::Integrate(Base range, T stepSize_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)) {}

    template<Scalar T>
    template<class Function>
    T Integrate<IntegrateMethod::Ladder, T, 1>::solve(Function func) const {
        const T& from = Base::from()[0];
        const T& to = Base::to()[0];
        T result = ((func(from) + func(to)) >> 1);
        T start(from + stepSize);
        while(start < to) {
            result += func(start);
            start += stepSize;
        }
        result *= stepSize;
        return result;
    }
    //////////////////////////////////Simpson//////////////////////////////////
    template<Scalar T>
    Integrate<IntegrateMethod::Simpson, T, 1>::Integrate(Base range, T stepSize_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)) {}

    template<Scalar T>
    template<class Function>
    T Integrate<IntegrateMethod::Simpson, T, 1>::solve(Function func) const {
        const T& from = Base::from()[0];
        const T& to = Base::to()[0];
        T result = func(from) + func(to);

        T odd = T(0);
        T even = T(0);
        bool b = true;
        T start = from + stepSize;
        while(start < to) {
            T& toChange = b ? odd : even;
            b = !b;
            toChange += func(start);
            start += stepSize;
        }
        odd <<= 2;
        even <<= 1;
        result += odd + even;
        result *= stepSize;
        result /= T(3);
        return result;
    }
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    /**
     * \param pointCount Points count on one side of the x-axis. There are same points on both sides.
     */
    template<Scalar T>
    Integrate<IntegrateMethod::Tanh_Sinh, T, 1>::Integrate(Base range, T stepSize_, uint64_t pointCount_)
            : Base(std::move(range)), stepSize(std::move(stepSize_)), pointCount(pointCount_) {}

    template<Scalar T>
    template<class Function>
    T Integrate<IntegrateMethod::Tanh_Sinh, T, 1>::solve(Function func) const {
        const T& from = Base::from()[0];
        const T& to = Base::to()[0];

        const T constant1 = (to - from) >> 1;
        const T constant2 = constant1 + from;
        const auto& PI_2 = T(M_PI_2);

        T result = PI_2 * func(constant2); //Integral value at t = 0
        T t = 0;
        for(uint64_t i = 0; i < pointCount; ++i) {
            t += stepSize;
            const T PI_2_sinh = PI_2 * sinh(t);
            const T cosh_PI_2_sinh = cosh(PI_2_sinh);
            const T phi = sinh(PI_2_sinh) / cosh_PI_2_sinh;
            const T phi_derivative = PI_2 * cosh(t) / square(cosh_PI_2_sinh);
            result += phi_derivative * (func(constant2 + constant1 * phi) + func(constant2 - constant1 * phi));
        }
        result *= constant1 * stepSize;
        return result;
    }
    //////////////////////////////////MonteCarlo//////////////////////////////////
    template<Scalar T, size_t dim>
    Integrate<IntegrateMethod::MonteCarlo, T, dim>::Integrate(Base range, uint64_t sampleCount_)
            : Base(std::move(range)), sampleCount(sampleCount_) {}

    template<Scalar T, size_t dim>
    template<class Function, class RandomGenerator>
    T Integrate<IntegrateMethod::MonteCarlo, T, dim>::solve(Function func, RandomGenerator& generator) const {
        T result = 0;
        for (uint64_t i = 0; i < sampleCount; ++i) {
            VectorType x = VectorType::template random_uniform(Base::from(), Base::to(), generator);
            toNextMean(result, i, func(x));
        }

        T factor = 1;
        for (size_t i = 0; i < Base::from().getLength(); ++i)
            factor *= Base::to()[i] - Base::from()[i];
        result *= factor;
        return result;
    }

    template<Scalar T, size_t dim>
    template<class Function, class RandomGenerator, class Executor>
    T Integrate<IntegrateMethod::MonteCarlo, T, dim>::parallel_solve(
            Function func,
            const Array<typename RandomGenerator::result_type>& seeds) const {
        const size_t numGen = seeds.getLength();
        Array<RandomGenerator> gens(numGen);
        for (size_t i = 0; i < numGen; ++i)
            gens[i] = RandomGenerator(seeds[i]);

        VectorND<T> results(numGen);
        auto future = Executor::parallel_for([this, func, &gens, &results](unsigned int i) { results[i] = solve(func, gens[i]); }, numGen, numGen);
        Executor::auto_wait(future);
        return mean(results);
    }

    template<Scalar T, size_t dim>
    template<class Function, class RandomGenerator>
    T Integrate<IntegrateMethod::MonteCarlo, T, dim>::solve_e(unsigned int numSequence, Function func, RandomGenerator& generator, T& deviation) const {
        assert(numSequence > 0);
        T mean = 0;
        T variance = 0;
        for (unsigned int i = 0; i < numSequence; ++i)
            toNextVariance(variance, mean, i, solve(func, generator));
        deviation = sqrt(variance);
        return mean;
    }

    template<Scalar T, size_t dim>
    template<class Functor1, class Functor2, class Distribution, class RandomGenerator>
    T Integrate<IntegrateMethod::MonteCarlo, T, dim>::solve(
            Functor1 func,
            Functor2 importance,
            Distribution& dist,
            RandomGenerator& generator) const {
        T result = 0;
        for (uint64_t i = 0; i < sampleCount; ++i) {
            const VectorType x = VectorType::template random_any(dim, dist, generator);
            toNextMean(result, i, func(x) / importance(x));
        }
        return result;
    }

    template<Scalar T, size_t dim>
    template<class Functor1, class Functor2, class Distribution, class RandomGenerator>
    T Integrate<IntegrateMethod::MonteCarlo, T, dim>::solve_e(
            unsigned int numSequence,
            Functor1 func,
            Functor2 importance,
            Distribution& dist,
            RandomGenerator& generator,
            T& deviation) const {
        assert(numSequence > 0);
        T mean = 0;
        T variance = 0;
        for (unsigned int i = 0; i < numSequence; ++i)
            toNextVariance(variance, mean, i, solve(func, importance, dist, generator));
        deviation = sqrt(variance);
        return mean;
    }
}
