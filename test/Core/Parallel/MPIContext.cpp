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
#include "Physica/Core/Parallel/MPIContext.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Test.h"

using namespace Physica;

namespace {
    void send() noexcept {
        const int id = MPIContext::getProcessID();
        Array<int, 1> arr{id};
        arr.send(0, 1);
        expect(arr[0] == 0);
    }

    void sendrecv() noexcept {
        const int id = MPIContext::getProcessID();
        int partner = id == 0 ? 1 : 0;
        Array<int, 1> arr{id};
        arr.sendrecv(partner, partner);
        expect(arr[0] == partner);
    }

    void bcast() noexcept {
        Array<int, 1> arr{MPIContext::getProcessID()};
        arr.bcast(0);
        expect(arr[0] == 0);
    }
}

int main() {
    std::ignore = MPIContext::getInstance();
    expect(MPIContext::getNumProcess() == 2);
    send();
    sendrecv();
    bcast();
    return 0;
}
