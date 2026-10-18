/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Utils/Container/ArrayND.h"
#include "Test.h"

using namespace Physica;

int main() {
    ArrayND<size_t, 2, 3, 4> arr{};
    arr.forND([](size_t& value, Index3D index) {
        value = index[0] * 100 + index[1] * 10 + index[2];
    });
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                expect(arr[i, j, k] == i * 100 + j * 10 + k);

    const auto& const_arr = arr;
    const_arr.forND([](const size_t& value, Index3D index) {
        expect(value == index[0] * 100 + index[1] * 10 + index[2]);
    });
    return 0;
}
