/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    void fiber() {
        auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        static_assert(x.fiber(1UZ, 2UZ, var()).isCompact());

        auto fiber = x.fiber(1UZ, var(), 2UZ);
        static_assert(!decltype(fiber)::isCompact());
        static_assert(decltype(fiber)::getStrideAtCompile() == Dynamic);
        expect(fiber.getStride() == x.getStride(1));
        for (size_t i = 0; i < fiber.getLength(); ++i)
            expect(x[1UZ, i, 2UZ] == fiber[i]);
    }

    void slice() {
        auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        auto slice = x.slice(1UZ, var(), var());
        static_assert(slice.isCompact());
        static_assert(slice.getMajor() == MatrixMajor::Row);
        for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
                expect(x[1UZ, r, c] == slice[r, c]);

        MatrixND<T> m = slice;
        expect(m == slice);
    }

    void block() {
        auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        const Index3D from{1, 0, 2};
        const Index3D count{2, 3, 1};
        auto b = x.block(from, count);
        static_assert(b.isStrided());

        expect(b.getShape() == count);
        expect(b.getStrides() == x.getStrides());
        expect(b.data_handle() == &x[1UZ, 0UZ, 2UZ]);
    }

    void flatten() {
        auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        auto flat = x.flatten();
        static_assert(flat.isCompact());

        expect(flat.getLength() == x.getSize());
        for (size_t i = 0; i < flat.getLength(); ++i)
            expect(flat[i] == x.calc(x.toIndexND(i)));

        flat = T(1);
        for (size_t i = 0; i < x.getSize(); ++i)
            expect(x.calc(x.toIndexND(i)) == T(1));
    }
}

int main() {
    fiber();
    slice();
    block();
    flatten();
    return 0;
}
