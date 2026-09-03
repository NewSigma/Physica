/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Parallel/MPI.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Test.h"

using namespace Physica;

namespace {
    void send() noexcept {
        const int rank = MPI::getRank();
        Array<int, 1> arr{rank};
        arr.pass(0, 1);
        expect(arr[0] == 0);
    }

    void sendrecv() noexcept {
        const int rank = MPI::getRank();
        int partner = rank == 0 ? 1 : 0;
        Array<int, 1> arr{rank};
        arr.sendrecv(partner, partner);
        expect(arr[0] == partner);
    }

    void bcast() noexcept {
        Array<int, 1> arr{MPI::getRank()};
        arr.bcast(0);
        expect(arr[0] == 0);
    }

    void reduce() noexcept {
        Array<int, 1> arr{MPI::getRank()};
        arr.reduce(0, MPI::ReduceOp::Sum);
        expect(arr[0] == 1);
    }

    void allreduce() noexcept {
        Array<int, 1> arr{MPI::getRank()};
        arr.allreduce(MPI::ReduceOp::Sum);
        expect(arr[0] == 1);
    }

    void gather() noexcept {
        const int rank = MPI::getRank();
        if (rank == 0) {
            Array<int, 2> arr{0, 0};
            arr.gather(0);
            expect(arr[0] == 0);
            expect(arr[1] == 1);
        }
        else {
            Array<int, 1> arr{rank};
            arr.gather(0);
        }
    }

    void scatter() noexcept {
        const int rank = MPI::getRank();
        if (rank == 0) {
            Array<int, 2> src{0, 1};
            src.scatter(0);
            expect(src[0] == 0);
            expect(src[1] == 1);
        }
        else {
            Array<int, 1> dst{100};
            dst.scatter(0);
            expect(dst[0] == 1);
        }
    }

    void allgather() noexcept {
        const int rank = MPI::getRank();
        Array<int, 2> arr{0, 0};
        arr[rank] = rank;
        arr.allgather();
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
