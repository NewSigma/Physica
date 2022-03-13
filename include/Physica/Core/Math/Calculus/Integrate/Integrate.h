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

#include "IntegrateRange.h"

namespace Physica::Core {
    enum IntegrateMethod {
        Rectangular,
        Ladder,
        Simpson,
        Tanh_Sinh
    };

    template<IntegrateMethod Method, class ScalarType, size_t dim>
    class Integrate;
    //////////////////////////////////Rectangular//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<Rectangular, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<Rectangular, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using Base::VectorType;
    private:
        ScalarType stepSize;
    public:
        explicit Integrate(Base range, ScalarType stepSize);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
    //////////////////////////////////Ladder//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<Ladder, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<Ladder, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using Base::VectorType;
    private:
        ScalarType stepSize;
    public:
        explicit Integrate(Base range, ScalarType stepSize);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
    //////////////////////////////////Simpson//////////////////////////////////
    template<class ScalarType, size_t dim>
    class Integrate<Simpson, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<Simpson, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using Base::VectorType;
    private:
        ScalarType stepSize;
    public:
        explicit Integrate(Base range, ScalarType stepSize);
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
    class Integrate<Tanh_Sinh, ScalarType, dim> : public IntegrateRange<ScalarType, dim> {
        using Base = IntegrateRange<ScalarType, dim>;
    };

    template<class ScalarType>
    class Integrate<Tanh_Sinh, ScalarType, 1> : public IntegrateRange<ScalarType, 1> {
        using Base = IntegrateRange<ScalarType, 1>;
    public:
        using Base::VectorType;
    private:
        ScalarType stepSize;
        uint64_t pointCount;
    public:
        Integrate(Base range, ScalarType stepSize, uint64_t pointCount);
        /* Operations */
        template<class Function>
        ScalarType solve(Function func) const;
    };
}

#include "IntegrateImpl.h"
