/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
