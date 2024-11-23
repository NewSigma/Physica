/*
 * Copyright 2020-2024 Weibo He.
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
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica::Core;
using RandomType = Random<MT19937, std::mt19937::default_seed>;

void crossProductTest() {
    using T = float32;
    VectorND<T> v1{3.845971,0.000000,0.000000};
    VectorND<T> v2{-0.007733,3.835502,0.000000};
    VectorND<T> v3(v1.crossProduct(v2));
    if (!scalarNear(v3.norm() / T(2), T(7.375614), 1E-7))
        exit(EXIT_FAILURE);
}

void hdfTest() {
#ifdef PHYSICA_HDF5
    /* Real */ {
        using T = float64;
        const auto data = VectorND<T>::template random_uniform<RandomType>(64);

        TempFile tmp("/tmp/tmpXXXXXX");
        {
            H5File h5f(tmp.getName(), H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat);
            data.write(h5f, "set");
        }
        VectorND<T> buffer(data.getLength());
        {
            H5File h5f(tmp.getName(), H5File::OpenFlag::ReadOnly);
            buffer.read(h5f, "set");
        }
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    /* Complex */ {
        using T = Complex<float64>;
        const auto data = VectorND<T>::template random_uniform<RandomType>(64);

        TempFile tmp("/tmp/tmpXXXXXX");
        {
            H5File h5f(tmp.getName(), H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat);
            data.write(h5f, "set");
        }
        VectorND<T> buffer(data.getLength());
        {
            H5File h5f(tmp.getName(), H5File::OpenFlag::ReadOnly);
            buffer.read(h5f, "set");
        }
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
#endif
}

int main() {
    crossProductTest();
    hdfTest();
    return 0;
}
