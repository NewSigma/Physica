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
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    /**
     * A continuous matrix is continuous in either row or column.
     */
    void continuousRowCol() noexcept {
        constexpr int Size = 8;
        using T = float32;
        auto& rng = RandomSource::getInstance();
        int r = std::uniform_int_distribution<>(0, Size - 1)(rng);
        int c = std::uniform_int_distribution<>(0, Size - 1)(rng);
        {
            using MatrixType = DenseMatrix<T, MatrixOption::Col>;
            const auto x = MatrixType(Size, Size);
            expect(x.data_ptr(r, c) == x.col(c).data() + r);
        }
        {
            using MatrixType = DenseMatrix<T, MatrixOption::Row>;
            const auto x = MatrixType(Size, Size);
            expect(x.data_ptr(r, c) == x.row(r).data() + c);
        }
    }

    void testHDF5() {
    #ifdef PHYSICA_HDF5
        /* Col float matrix */ {
            using T = float32;
            using MatrixType = DenseMatrix<T>;
            const auto data = MatrixType::random_uniform<RandomSource>(32, 16);

            TempFile tmp("/tmp/tmpXXXXXX");
            auto h5f = H5File::open(tmp.getName());
            data.write(h5f, "/set");

            MatrixType buffer(data.getRow(), data.getCol());
            buffer.read(h5f, "/set");
            expect(data == buffer);
        }
        /* Row double matrix */ {
            using T = float64;
            using MatrixType = DenseMatrix<T, MatrixOption::Row>;
            const auto data = MatrixType::random_uniform<RandomSource>(16, 20);

            TempFile tmp("/tmp/tmpXXXXXX");
            auto h5f = H5File::open(tmp.getName());
            data.write(h5f, "/set");

            MatrixType buffer(data.getRow(), data.getCol());
            buffer.read(h5f, "/set");
            expect(data == buffer);
        }
        /* Row complex float matrix */ {
            using T = Complex<float32>;
            using MatrixType = DenseMatrix<T, MatrixOption::Row>;
            const auto data = MatrixType::random_uniform<RandomSource>(16, 12);

            TempFile tmp("/tmp/tmpXXXXXX");
            auto h5f = H5File::open(tmp.getName());
            data.write(h5f, "/set");

            MatrixType buffer(data.getRow(), data.getCol());
            buffer.read(h5f, "/set");
            expect(data == buffer);
        }
    #endif
    }
}

int main() {
    continuousRowCol();
    testHDF5();
    return 0;
}
