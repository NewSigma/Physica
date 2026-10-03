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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void struct_bind() noexcept {
        Vector3D<Diff<float64, DiffMode::Forward>> arr = Vector3D<float64>{1, 2, 3};
        auto [x, y, z] = arr;
        expect(x == 1);
        expect(y == 2);
        expect(z == 3);
    }

    void forward() {
        syntax_only([]() {
            // Test mixed types for binary operation
            using V0 = VectorND<float64>;
            using V1 = VectorND<Diff<float64, DiffMode::Forward>>;
            auto v0 = V0::random_uniform<RandomSource>(16);
            auto v1 = V1::random_uniform<RandomSource>(16);
            V1 xx = hadamard(v0, v1);
        });
    }

    void reverse() {
        using VectorType = VectorND<Diff<float64, DiffMode::Reverse, 1>>;
        auto v = VectorType::random_uniform<RandomSource>(16);
        v.sum().reverse();
        for (const auto& elem : v)
            expect(elem.grad() == float64(1));
    }

    void reverse_sum_unary() {
        using V = VectorND<Diff<float64, DiffMode::Reverse, 1>>;
        auto v = V::random_uniform<RandomSource>(16);
        const auto check = [&](const auto& func) {
            v.zero_grad();
            func(v).sum().reverse();
            for (const auto& x : v) {
                using dfloat = Diff<float64, DiffMode::Forward, 1>;
                const dfloat ref = func(dfloat(x.value(), float64(1)));
                expect<RandomSource>(scalarNear(x.grad(), ref.grad(), 1E-8));
            }
        };

        check([](const auto& x) { return exp(x); });
        check([](const auto& x) { return expm1(x); });
        check([](const auto& x) { return sqrt(x); });
        check([](const auto& x) { return cbrt(x); });
        check([](const auto& x) { return sin(x); });
        check([](const auto& x) { return cos(x); });
        check([](const auto& x) { return tan(x); });
        check([](const auto& x) { return sec(x); });
        check([](const auto& x) { return tanh(x); });
        check([](const auto& x) { return cosh(x); });
        check([](const auto& x) { return sech(x); });
        check([](const auto& x) { return sigmoid(x); });
        check([](const auto& x) { return softplus(x); });
        check([](const auto& x) { return ln(x); });
        check([](const auto& x) { return ln1p(x); });
        check([](const auto& x) { return lncosh(x); });
        check([](const auto& x) { return arcsinh(x); });
        check([](const auto& x) { return arctanh(x); });
        check([](const auto& x) { return reciprocal(x); });
        check([](const auto& x) { return relu(x); });
        check([](const auto& x) { return abs(x); });
        check([](const auto& x) { return square(x); });
        check([](const auto& x) { return pow(x, float64(3)); });
    }

    void reverse_softmax_sum() {
        using VectorType = VectorND<Diff<float64, DiffMode::Reverse, 1>>;
        auto v = VectorType::random_uniform<RandomSource>(16);
        float64 sum = softmax(v).sum().reverse();
        expect<RandomSource>(scalarNear(sum, float64(1), 1E-8));
        for (const auto& elem : v)
            expect<RandomSource>(scalarNear(elem.grad(), float64(0), 1E-8));
    }

    void reverse_mean() {
        const float64 factor = float64(1) / 16;
        auto v = VectorND<Diff<float64, DiffMode::Reverse, 1>>::random_uniform<RandomSource>(16);
        v.mean().reverse();
        for (const auto& elem : v)
            expect<RandomSource>(scalarNear(elem.grad(), factor, 1E-8));

        v.zero_grad();
        v.mean_stable().reverse();
        for (const auto& elem : v)
            expect<RandomSource>(scalarNear(elem.grad(), factor, 1E-8));
    }

    void reverse_variance() {
        auto v = VectorND<Diff<float64, DiffMode::Reverse, 1>>::random_uniform<RandomSource>(16);
        v.variance().reverse();
        const float64 mean = v.values().mean();
        for (const auto& elem : v) {
            const float64 x = elem.value().value();
            expect<RandomSource>(scalarNear(elem.grad(), float64(2) * (x - mean) / 16, 1E-8));
        }
    }

    void reverse_pow_vector() {
        // Covers grad(vector) path for pow()
        auto x = VectorND<Diff<float64, DiffMode::Reverse, 1>>::random_uniform<RandomSource>(16);
        const auto w = VectorND<float64>::random_uniform<RandomSource>(16);
        const auto expr = pow(x, float64(3));
        expr.reverse(w);
        for (auto [xi, wi] : zip(x, w)) {
            const float64 value = xi.value().value();
            expect<RandomSource>(scalarNear(xi.grad(), float64(3) * wi * square(value), 1E-8));
        }
    }

    void test_hdf5() {
    #ifdef PHYSICA_HDF5
        using dfloat = Diff<float64, DiffMode::Forward, 1>;
        auto data = VectorND<dfloat>::random_uniform<RandomSource>(36);
        data.grads().random_uniform<RandomSource>();

        TempFile tmp("/tmp/tmpXXXXXX");
        {
            auto h5f = H5File::open(tmp.getName());
            data.write(h5f, "x");
        }
        VectorND<dfloat> result;
        auto h5f = H5File::open(tmp.getName(), H5File::ReadOnly);
        result.read(h5f, "x");
        expect(data == result);
    #endif
    }
}

int main() {
    struct_bind();
    forward();
    reverse();
    reverse_sum_unary();
    reverse_softmax_sum();
    reverse_mean();
    reverse_variance();
    reverse_pow_vector();
    test_hdf5();
    return 0;
}
