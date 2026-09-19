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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tensor4x4x4 = DenseTensor<T, 4, 4, 4>;
using RandomSource = Random<>;

namespace {
    void compact() {
        const Tensor3D<T> x = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        expect(x.data() == x.asArray().data());
        expect(x.data() == x.data_handle());
        expect(x.data_ptr({1, 2, 3}) == x.data() + x.toIndex1D({1, 2, 3}));

        expect(x.getStride(0) == 16);
        expect(x.getStride(1) == 4);
        expect(x.getStride(2) == 1);

        for (size_t i = 0; i < x.dim(0); ++i)
            for (size_t j = 0; j < x.dim(1); ++j)
                for (size_t k = 0; k < x.dim(2); ++k)
                    expect(x.data_ptr({i, j, k}) == &x[i, j, k]);
    }
}

static_assert(Tensor3D<T>::isCompact(), "DenseTensor is a compact object");
static_assert(Tensor3D<T>::isStrided(), "DenseTensor is a strided object");
static_assert(Tensor3D<T>::getStrideAtCompile()[0] == Dynamic, "Dynamic shape strides are unknown at compile time");
static_assert(Tensor3D<T>::getStrideAtCompile()[2] == 1, "DenseTensor is row-major");
static_assert(Tensor4x4x4::getStrideAtCompile()[0] == 16);
static_assert(Tensor4x4x4::getStrideAtCompile()[1] == 4);
static_assert(Tensor4x4x4::getStrideAtCompile()[2] == 1);

static_assert(Tensor3D<T>::getStrideAtCompile(0) == Tensor3D<T>::getStrideAtCompile()[0]);
static_assert(Tensor3D<T>::getStrideAtCompile(2) == Tensor3D<T>::getStrideAtCompile()[2]);
static_assert(Tensor4x4x4::getStrideAtCompile(0) == Tensor4x4x4::getStrideAtCompile()[0]);
static_assert(Tensor4x4x4::getStrideAtCompile(1) == Tensor4x4x4::getStrideAtCompile()[1]);
static_assert(Tensor4x4x4::getStrideAtCompile(2) == Tensor4x4x4::getStrideAtCompile()[2]);

int main() {
    compact();
    return 0;
}
