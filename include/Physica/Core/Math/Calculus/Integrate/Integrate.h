/*
 * Copyright 2020-2025 Weibo He.
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

#include "IntegrateImpl/IntegrateRange.h"

namespace Physica {
    enum class IntegrateMethod {
        Rectangular,
        Ladder,
        Simpson,
        Tanh_Sinh,
        MonteCarlo
    };

    template<IntegrateMethod Method, Scalar T, size_t dim>
    class Integrate;
    //////////////////////////////////Rectangular//////////////////////////////////
    template<Scalar T, size_t dim>
    class Integrate<IntegrateMethod::Rectangular, T, dim> : public IntegrateRange<T, dim> {
        using Base = IntegrateRange<T, dim>;
    };

    template<Scalar T>
    class Integrate<IntegrateMethod::Rectangular, T, 1> : public IntegrateRange<T, 1> {
        using Base = IntegrateRange<T, 1>;
    public:
        using typename Base::VectorType;
    private:
        T stepSize;
    public:
        explicit Integrate(Base range, T stepSize_);
        /* Operations */
        template<class Function>
        T solve(Function func) const;
    };
    //////////////////////////////////Ladder//////////////////////////////////
    template<Scalar T, size_t dim>
    class Integrate<IntegrateMethod::Ladder, T, dim> : public IntegrateRange<T, dim> {
        using Base = IntegrateRange<T, dim>;
    };

    template<Scalar T>
    class Integrate<IntegrateMethod::Ladder, T, 1> : public IntegrateRange<T, 1> {
        using Base = IntegrateRange<T, 1>;
    public:
        using typename Base::VectorType;
    private:
        T stepSize;
    public:
        explicit Integrate(Base range, T stepSize_);
        /* Operations */
        template<class Function>
        T solve(Function func) const;
    };
    //////////////////////////////////Simpson//////////////////////////////////
    template<Scalar T, size_t dim>
    class Integrate<IntegrateMethod::Simpson, T, dim> : public IntegrateRange<T, dim> {
        using Base = IntegrateRange<T, dim>;
    };

    template<Scalar T>
    class Integrate<IntegrateMethod::Simpson, T, 1> : public IntegrateRange<T, 1> {
        using Base = IntegrateRange<T, 1>;
    public:
        using typename Base::VectorType;
    private:
        T stepSize;
    public:
        explicit Integrate(Base range, T stepSize_);
        /* Operations */
        template<class Function>
        T solve(Function func) const;
    };
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    /**
     * Reference:
     * [1] arXiv:2007.15057; https://doi.org/10.48550/arXiv.2007.15057
     */
    template<Scalar T, size_t dim>
    class Integrate<IntegrateMethod::Tanh_Sinh, T, dim> : public IntegrateRange<T, dim> {
        using Base = IntegrateRange<T, dim>;
    };

    template<Scalar T>
    class Integrate<IntegrateMethod::Tanh_Sinh, T, 1> : public IntegrateRange<T, 1> {
        using Base = IntegrateRange<T, 1>;
    public:
        using typename Base::VectorType;
    private:
        T stepSize;
        uint64_t pointCount;
    public:
        Integrate(Base range, T stepSize_, uint64_t pointCount_);
        /* Operations */
        template<class Function>
        T solve(Function func) const;
    };
    //////////////////////////////////MonteCarlo//////////////////////////////////
    template<Scalar T, size_t dim>
    class Integrate<IntegrateMethod::MonteCarlo, T, dim> : public IntegrateRange<T, dim> {
        using Base = IntegrateRange<T, dim>;
    public:
        using typename Base::VectorType;
    private:
        uint64_t sampleCount;
    public:
        Integrate(Base range, uint64_t sampleCount_);
        /* Operations */
        template<class Function, RNG R>
        T solve(Function func) const;
        template<class Function, RNG R>
        T solve_e(unsigned int numSequence, Function func, T& deviation) const;

        template<class Functor1, class Functor2, class Distribution, RNG R>
        T solve(Functor1 func, Functor2 importance, Distribution& disterator) const;
        template<class Functor1, class Functor2, class Distribution, RNG R>
        T solve_e(unsigned int numSequence,
                           Functor1 func,
                           Functor2 importance,
                           Distribution& dist,
                           T& deviation) const;
    };
}

#include "IntegrateImpl/IntegrateImpl.h"
