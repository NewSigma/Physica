/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "IntegrateImpl/IntegrateRange.h"

namespace Physica {
    enum class IntegrateMethod : int8_t {
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
        [[nodiscard]] T solve(std::invocable<T> auto fn) const;
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
        [[nodiscard]] T solve(std::invocable<T> auto fn) const;
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
        [[nodiscard]] T solve(std::invocable<T> auto fn) const;
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
        [[nodiscard]] T solve(std::invocable<T> auto fn) const;
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
        template<RNG R>
        [[nodiscard]] T solve(std::invocable<VectorType> auto fn) const;
        template<RNG R>
        [[nodiscard]] T solve_e(unsigned int numSequence, std::invocable<VectorType> auto fn, T& deviation) const;

        template<RNG R>
        [[nodiscard]] T solve(std::invocable<VectorType> auto fn, std::invocable<VectorType> auto importance, auto& distribution) const;
        template<RNG R>
        [[nodiscard]] T solve_e(unsigned int numSequence,
                                std::invocable<VectorType> auto fn,
                                std::invocable<VectorType> auto importance,
                                auto& distribution,
                                T& deviation) const;
    };
}

#include "IntegrateImpl/IntegrateImpl.h"
