/*
 * Copyright 2020-2022 WeiBo He.
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
#include <fstream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Utils;

void crossProductTest() {
    using T = Scalar<Float>;
    Vector<T> v1{3.845971,0.000000,0.000000};
    Vector<T> v2{-0.007733,3.835502,0.000000};
    Vector<T> v3(v1.crossProduct(v2));
    if (!scalarNear(v3.norm() / T(2), T(7.375614), 1E-7))
        exit(EXIT_FAILURE);
}

void hdfTest() {
    std::mt19937 gen{};
    /* Real */ {
        using T = Scalar<Double>;
        const auto data = Vector<T>::random_uniform(64, gen);

        TempFile tmp("/tmp/tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        auto space = H5DataSpace<1>::makeDataSpace({data.getLength()});
        auto dataset = h5f.createDataSet("/set", T::getH5DataType(), space);
        data.write(dataset, space);

        Vector<T> buffer(data.getLength());
        buffer.read(dataset, space);
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    /* Complex */ {
        using T = ComplexScalar<Scalar<Double>>;
        const auto data = Vector<T>::random_uniform(64, gen);

        TempFile tmp("/tmp/tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        auto space = H5DataSpace<1>::makeDataSpace({data.getLength()});
        auto dataset = h5f.createDataSet("/set", T::getH5DataType(), space);
        data.write(dataset, space);

        Vector<T> buffer(data.getLength());
        buffer.read(dataset, space);
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
}

int main() {
    crossProductTest();
    hdfTest();
    return 0;
}
