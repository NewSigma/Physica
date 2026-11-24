/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <chrono>
#include <thread>
#include "Physica/Core/Parallel/Algorithm/Sequential.h"
#include "Physica/Core/Parallel/Algorithm/Thread.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using MatrixType = DenseMatrix<ScalarType>;
using VectorType = VectorND<ScalarType>;
using RandomSource = Random<>;

namespace {
    void parallel_for_test() {
        const MatrixType A = MatrixType::random_uniform<RandomSource>(4);
        ThreadPool::numThreadRequired = 4;

        VectorType result_seq(4);
        parallel_for<Sequential>([&](unsigned int i) {
            result_seq[i] = A.col(i).sum();
        }, 4, 4);

        VectorType result_par(4);
        parallel_for<Thread>([&](unsigned int i) {
            result_par[i] = A.col(i).sum();
        }, 4, 4).wait();

        expect(result_seq == result_par);
    }

    void resize_schedule() {
        // Regression test that using noop_coroutine as the sentinel does not crash
        ThreadPool& pool = ThreadPool::getInstance();
        pool.resize(2);
        expect(pool.getNumThreads() == 2);
        []() -> Task {
            // optnone is required to reproduce the crash
            auto schedule_optnone = [][[clang::optnone]]() -> Task {
                co_await ThreadPool::getInstance();
                // Sleep until the parent reaches await_ready()
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            };

            co_await schedule_optnone();
        }().wait();
        pool.resize(4);
        expect(pool.getNumThreads() == 4);
    }
}

int main() {
    parallel_for_test();
    resize_schedule();
    return 0;
}
