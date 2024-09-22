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

#include "IntegrateRange.h"

namespace Physica::Core {
    enum class IntegrateMethod {
        Rectangular,
        Ladder,
        Simpson,
        Tanh_Sinh,
        MonteCarlo
    };

    template<IntegrateMethod Method, class ScalarType, size_t dim>
    class Integrate;
    //////////////////////////////////Rectangular//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<IntegrateMethod::Rectangular, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<IntegrateMethod::Rectangular, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using typename Base::VectorType;
    private:
        ScalarType stepSize;
    public:
        explicit Integrate(Base range, ScalarType stepSize_);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
    //////////////////////////////////Ladder//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<IntegrateMethod::Ladder, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<IntegrateMethod::Ladder, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using typename Base::VectorType;
    private:
        ScalarType stepSize;
    public:
        explicit Integrate(Base range, ScalarType stepSize_);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
    //////////////////////////////////Simpson//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<IntegrateMethod::Simpson, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<IntegrateMethod::Simpson, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using typename Base::VectorType;
    private:
        ScalarType stepSize;
    public:
        explicit Integrate(Base range, ScalarType stepSize_);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    /**
     * Reference:
     * [1] Vanherck, Joren Sorée, Bart Magnus, Wim.
     * Tanh-sinh quadrature for single and multiple integration using floating-point arithmetic. arXiv:2007.15057
     */
    template<class ScalarType, size_t dim>
    class Integrate<IntegrateMethod::Tanh_Sinh, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<IntegrateMethod::Tanh_Sinh, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using typename Base::VectorType;
    private:
        ScalarType stepSize;
        uint64_t pointCount;
    public:
        Integrate(Base range, ScalarType stepSize_, uint64_t pointCount_);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
    //////////////////////////////////MonteCarlo//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<IntegrateMethod::MonteCarlo, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    public:
        using typename Base::VectorType;
    private:
        uint64_t sampleCount;
    public:
        Integrate(Base range, uint64_t sampleCount_);
        /* Operations */
        template<class Function, class RandomGenerator>
        ScalarType solve(Function func, RandomGenerator& generator) const;
        template<class Function, class RandomGenerator, class Executor>
        ScalarType parallel_solve(Function func, const Array<typename RandomGenerator::result_type>& seeds) const;
        template<class Function, class RandomGenerator>
        ScalarType solve_e(unsigned int numSequence, Function func, RandomGenerator& generator, ScalarType& deviation) const;

        template<class Functor1, class Functor2, class Distribution, class RandomGenerator>
        ScalarType solve(Functor1 func, Functor2 importance, Distribution& dist, RandomGenerator& generator) const;
        template<class Functor1, class Functor2, class Distribution, class RandomGenerator>
        ScalarType solve_e(unsigned int numSequence,
                           Functor1 func,
                           Functor2 importance,
                           Distribution& dist,
                           RandomGenerator& generator,
                           ScalarType& deviation) const;
    };
}

#include "IntegrateImpl.h"
