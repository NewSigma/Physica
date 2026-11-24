/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Parallel/MPI.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Test.h"

using namespace Physica;

namespace {
    void send() noexcept {
        const int rank = MPI::getRank();
        Array<int, 1> arr{rank};
        MPI::pass(0, 1, arr).wait();
        expect(arr[0] == 0);
    }

    void sendrecv() noexcept {
        const int rank = MPI::getRank();
        int partner = rank == 0 ? 1 : 0;
        Array<int, 1> arr{rank};
        MPI::sendrecv(partner, partner, arr);
        expect(arr[0] == partner);
    }

    void bcast() noexcept {
        Array<int, 1> arr{MPI::getRank()};
        MPI::bcast(0, arr).wait();
        expect(arr[0] == 0);
    }

    void reduce() noexcept {
        Array<int, 1> arr{MPI::getRank()};
        MPI::reduce(0, arr, arr, MPI::ReduceOp::Sum).wait();
        expect(arr[0] == 1);
    }

    void allreduce() noexcept {
        Array<int, 1> arr{MPI::getRank()};
        MPI::allreduce(arr, arr, MPI::ReduceOp::Sum).wait();
        expect(arr[0] == 1);
    }

    void gather() noexcept {
        const int rank = MPI::getRank();
        if (rank == 0) {
            Array<int, 2> arr{0, 0};
            MPI::gather(0, arr, arr).wait();
            expect(arr[0] == 0);
            expect(arr[1] == 1);
        }
        else {
            Array<int, 1> arr{rank};
            MPI::gather(0, arr, arr).wait();
        }
    }

    void scatter() noexcept {
        const int rank = MPI::getRank();
        if (rank == 0) {
            Array<int, 2> src{0, 1};
            MPI::scatter(0, src, src).wait();
            expect(src[0] == 0);
            expect(src[1] == 1);
        }
        else {
            Array<int, 1> dst{100};
            MPI::scatter(0, dst, dst).wait();
            expect(dst[0] == 1);
        }
    }

    void allgather() noexcept {
        const int rank = MPI::getRank();
        Array<int, 2> arr{0, 0};
        arr[rank] = rank;
        MPI::allgather(arr, arr).wait();
        expect(arr[0] == 0);
        expect(arr[1] == 1);
    }
}

int main() {
    std::ignore = MPI::getInstance();
    expect(MPI::getNumRank() == 2);
    send();
    sendrecv();
    bcast();
    reduce();
    allreduce();
    gather();
    scatter();
    allgather();
    return 0;
}
