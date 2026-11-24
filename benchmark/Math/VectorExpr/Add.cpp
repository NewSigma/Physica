/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

// Benchmark that we can recognize fma patterns
namespace {
    constexpr size_t size = 1024;

    template<Scalar T>
    void add_scalar1(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> b = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = hadamard(a, b) + T(0.1);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    template<Scalar T>
    void add_scalar2(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> c = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = a * T(0.1) + c;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    template<Scalar T>
    void add_scalar3(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = a * T(0.1) + T(0.2);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    template<Scalar T>
    void add_scalar_inplace(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = a * T(0.1);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_add(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    template<Scalar T>
    void add_vector(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> b = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> c = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = hadamard(a, b) + c;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    template<Scalar T>
    void add_vector_inplace(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> b = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = hadamard(a, b);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_add(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(add_scalar1<float64>)->Name("add scalar1");
BENCHMARK(add_scalar2<float64>)->Name("add scalar2");
BENCHMARK(add_scalar3<float64>)->Name("add scalar3");
BENCHMARK(add_scalar_inplace<float64>)->Name("add scalar inplace");
BENCHMARK(add_vector<float64>)->Name("add vector");
BENCHMARK(add_vector_inplace<float64>)->Name("add vector inplace");
