/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Random/Sobol.h"

using namespace Physica;

Sobol::Sobol() : numStep(0), mask(MaxDim * MaxBit), buffer(MaxDim, 0) {
    constexpr static int MaskInit[]{1, 1, 1, 1, 1, 1, 3, 1, 3, 3, 1, 1, 5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9};
    for (unsigned int i = 0; i < sizeof(MaskInit) / sizeof(int); ++i)
        mask[i] = MaskInit[i];

    Array<unsigned int*> iu(MaxBit);
    for (int j = 0, k = 0; j < MaxBit; ++j, k += MaxDim)
        iu[j] = &mask[k];

    for (int k = 0; k < MaxDim; k++) {
        for (int j = 0; j < Degree[k]; ++j)
            iu[j][k] <<= (MaxBit - 1 - j);

        int p = PolyCoeff[k];
        for (int j = Degree[k]; j < MaxBit; ++j) {
            unsigned int i = iu[j - Degree[k]][k];
            i ^= (i >> Degree[k]);
            for (int l = Degree[k] - 1; l >= 1; --l) {
                if (p % 2 == 1)
                    i ^= iu[j - l][k];
                p >>= 1;
            }
            iu[j][k] = i;
        }
    }
}

void Sobol::step(int i) {
    for (int _ = 0; _ < i; ++_) {
        const int m = pre_step();
        for (int j = 0; j < MaxDim; ++j)
            buffer[j] ^= mask[m + j];
    }
}

void Sobol::reset() {
    numStep = 0;
    for (auto& x : buffer)
        x = 0;
}

void Sobol::swap(This& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(numStep, obj.numStep);
    mask.swap(obj.mask);
    buffer.swap(obj.buffer);
}

int Sobol::pre_step() {
    unsigned int im = numStep++;
    int j = 0;
    for (; j < MaxBit; ++j) {
        if (im % 2 == 0)
            break;
        im >>= 1;
    }
    assert(j < MaxBit && "MaxBit too small in sobol");
    return j * MaxDim;
}

auto Sobol::getInstance() noexcept -> This& {
    thread_local static This instance{};
    return instance;
}
