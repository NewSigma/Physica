/*
 * Copyright 2024 Weibo He.
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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Utils/Unix/TempFile.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

static void test_rev() {
    using VectorType = VectorND<Diff<float64, DiffMode::Reverse, 1>>;
    auto v = VectorType::random_uniform<RandomSource>(16);
    v.sum().reverse();
    for (size_t i = 0; i < v.getLength(); ++i)
        if (v[i].grad() != float64(1)) [[unlikely]]
            exit(EXIT_FAILURE);
}

static void test_hdf5() {
#ifdef PHYSICA_HDF5
    using dfloat = Diff<float64, DiffMode::Forward, 1>;
    auto data = VectorND<dfloat>::random_uniform<RandomSource>(36);
    data.grads().random_uniform<RandomSource>();

    TempFile tmp("/tmp/tmpXXXXXX");
    {
        auto h5f = H5File::create(tmp.getName());
        data.write(h5f, "x");
    }
    VectorND<dfloat> result;
    auto h5f = H5File(tmp.getName());
    result.read(h5f, "x");
    if (data != result)
        exit(EXIT_FAILURE);
#endif
}

int main() {
    test_rev();
    test_hdf5();
    return 0;
}
