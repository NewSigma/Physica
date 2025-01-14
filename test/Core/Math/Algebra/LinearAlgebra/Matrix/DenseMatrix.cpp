/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Utils/Unix/TempFile.h"

using namespace Physica;
using RandomType = Random<MT19937>;

int main() {
#ifdef PHYSICA_HDF5
    /* Col float matrix */ {
        using T = float32;
        using MatrixType = DenseMatrix<T>;
        const auto data = MatrixType::random_uniform<RandomType>(32, 16);

        TempFile tmp("/tmp/tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        data.write(h5f, "/set");

        MatrixType buffer(data.getRow(), data.getCol());
        buffer.read(h5f, "/set");
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    /* Row double matrix */ {
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Row | MatrixOption::Vector>;
        const auto data = MatrixType::random_uniform<RandomType>(16, 20);

        TempFile tmp("/tmp/tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        data.write(h5f, "/set");

        MatrixType buffer(data.getRow(), data.getCol());
        buffer.read(h5f, "/set");
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
    /* Row complex float matrix */ {
        using T = Complex<float32>;
        using MatrixType = DenseMatrix<T, MatrixOption::Row | MatrixOption::Vector>;
        const auto data = MatrixType::random_uniform<RandomType>(16, 12);

        TempFile tmp("/tmp/tmpXXXXXX");
        H5File h5f(tmp.getName(), H5File::OpenFlag(H5File::OpenFlag::ReadWrite | H5File::OpenFlag::Creat));
        data.write(h5f, "/set");

        MatrixType buffer(data.getRow(), data.getCol());
        buffer.read(h5f, "/set");
        if (data != buffer)
            exit(EXIT_FAILURE);
    }
#endif
    return 0;
}
