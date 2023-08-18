/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Utils;

int main() {
    std::mt19937 gen{};
    /* Column float matrix */ {
        using T = Scalar<Float>;
        using MatrixType = DenseMatrix<T>;
        const auto data = MatrixType::random_uniform(32, 16, gen);

        TempFile tmp("tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        auto space = H5DataSpace<2>::makeDataSpace({data.getRow(), data.getColumn()});
        auto dataset = h5f.createDataSet("/set", T::getH5DataType(), space);
        data.write(dataset, space);

        MatrixType buffer(data.getRow(), data.getColumn());
        buffer.read(dataset, space);
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    /* Row double matrix */ {
        using T = Scalar<Double>;
        using MatrixType = DenseMatrix<T, MatrixOption::Row | MatrixOption::Vector>;
        const auto data = MatrixType::random_uniform(16, 20, gen);

        TempFile tmp("tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        auto space = H5DataSpace<2>::makeDataSpace({data.getRow(), data.getColumn()});
        auto dataset = h5f.createDataSet("/set", T::getH5DataType(), space);
        data.write(dataset, space);

        MatrixType buffer(data.getRow(), data.getColumn());
        buffer.read(dataset, space);
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    /* Row complex float matrix */ {
        using T = ComplexScalar<Scalar<Float>>;
        using MatrixType = DenseMatrix<T, MatrixOption::Row | MatrixOption::Vector>;
        const auto data = MatrixType::random_uniform(16, 12, gen);

        TempFile tmp("tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        auto space = H5DataSpace<2>::makeDataSpace({data.getRow(), data.getColumn()});
        auto dataset = h5f.createDataSet("/set", T::getH5DataType(), space);
        data.write(dataset, space);

        MatrixType buffer(data.getRow(), data.getColumn());
        buffer.read(dataset, space);
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    return 0;
}
