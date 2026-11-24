/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Test.h"

using namespace Physica;

namespace {
    void fiber(const DenseTensor<float64, 3>& x) {
        auto fiber = x.fiber(1UZ, var(), 2UZ);
        for (size_t i = 0; i < fiber.getLength(); ++i)
            expect(x[1UZ, i, 2UZ] == fiber[i]);

        VectorND<float64> v = fiber;
        expect(v == fiber);
    }

    void slice(const DenseTensor<float64, 3>& x) {
        auto slice = x.slice(1UZ, var(), var());
        for (size_t r = 0; r < x.dim(1); ++r)
            for (size_t c = 0; c < x.dim(2); ++c)
                expect(x[1UZ, r, c] == slice[r, c]);

        MatrixND<float64> m = slice;
        expect(m == slice);
    }

    template<Tensor X>
    void block(const X& x, typename X::IndexType from, typename X::IndexType count) {
        auto b = x.block(from, count);
        expect(b.getShape() == count);
        for (int i = 0; i < X::ndim(); ++i)
            expect(b.dim(i) == count[i]);

        for (size_t i = 0; i < b.getSize(); ++i) {
            const auto local = b.toIndexND(i);
            auto global = local;
            for (auto&& [g, f] : zip(global, from))
                g += f;
            expect(x[global] == b[local]);
        }

        X copy = b;
        for (size_t i = 0; i < b.getSize(); ++i) {
            const auto idx = b.toIndexND(i);
            expect(copy[idx] == b[idx]);
        }
    }
}

int main() {
    auto x = DenseTensor<float64, 3>::random_uniform<Random<>>({4, 4, 4});
    fiber(x);
    slice(x);
    block(x, Index3D{1, 0, 2}, Index3D{2, 3, 1});

    auto x4 = DenseTensor<float64, 4>::random_uniform<Random<>>({4, 2, 3, 2});
    block(x4, Index4D{1, 0, 2, 0}, Index4D{2, 1, 1, 2});
    return 0;
}
