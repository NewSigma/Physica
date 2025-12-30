/*
 * Copyright 2020-2025 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

namespace {
    template<Vector V>
    consteval void rangeTest() noexcept {
        using It = PtrIteratorF<V>;
        static_assert(std::indirectly_readable<It>);
        static_assert(std::indirectly_writable<It, long>);
        static_assert(std::incrementable<It>);
        static_assert(std::sized_sentinel_for<It, It>);
        static_assert(std::contiguous_iterator<It>);

        static_assert(std::ranges::sized_range<V>);
        static_assert(std::ranges::contiguous_range<V>);
        static_assert(std::ranges::common_range<V>);
        static_assert(std::ranges::viewable_range<V>);
    }

    void strucBindTest() noexcept {
        Vector3D<float64> arr{1, 2, 3};
        auto [x, y, z] = arr;
        expect(x == 1);
        expect(y == 2);
        expect(z == 3);
    }

    void crossProductTest() {
        using T = float32;
        VectorND<T> v1{3.845971, 0.000000, 0.000000};
        VectorND<T> v2{-0.007733, 3.835502, 0.000000};
        VectorND<T> v3(v1.cross(v2));
        expect(scalarNear(v3.norm() / T(2), T(7.375614), 1E-7));
    }

    void innerDotTest() {
        using T = float32;
        const auto v1 = VectorND<Complex<T>>::random_uniform<RandomSource>(16);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(16);
        const auto dot1 = v1 * v2;
        const auto dot2 = hadamard(v1, v2).sum();
        expect(scalarNear(dot1, dot2, 1E-6));
    }

    void hdfTest() {
    #ifdef PHYSICA_HDF5
        /* Real */ {
            using T = float64;
            const auto data = VectorND<T>::template random_uniform<RandomSource>(64);

            TempFile tmp("/tmp/tmpXXXXXX");
            {
                auto h5f = H5File::open(tmp.getName());
                data.write(h5f, "set");
            }
            VectorND<T> buffer(data.getLength());
            {
                auto h5f = H5File::open(tmp.getName(), H5File::OpenFlag::ReadOnly);
                buffer.read(h5f, "set");
            }
            expect(data == buffer);
        }
        /* Complex */ {
            using T = Complex<float64>;
            const auto data = VectorND<T>::template random_uniform<RandomSource>(64);

            TempFile tmp("/tmp/tmpXXXXXX");
            {
                auto h5f = H5File::open(tmp.getName());
                data.write(h5f, "set");
            }
            VectorND<T> buffer(data.getLength());
            {
                auto h5f = H5File::open(tmp.getName(), H5File::OpenFlag::ReadOnly);
                buffer.read(h5f, "set");
            }
            expect(data == buffer);
        }
    #endif
    }

    void lnSumExpTest() {
        /* Complex overflow test */ {
            VectorND<cfloat64> v{-1071, -739};
            expect(v.lnSumExp().isFinite());
        }
        /* Diff test */ {
            using dfloat = Diff<float32, DiffMode::Reverse, 1>;
            const auto x = VectorND<dfloat>::random_uniform<RandomSource>(8);
            x.lnSumExp().reverse();

            for (size_t i = 0; i < x.getLength(); ++i)
                expect(scalarNear(x.grads()[i], x.values().softmax(i), 1E-6));
        }
    }

    void crossEntropyTest() {
        /* Select test */ {
            const VectorND<float32> result{-3.34036088, -109.5531235, 13.51656151, 11.29175949};
            const float32 l1 = result.crossEntropy(3);
            expect(l1.isFinite());

            const float32 l2 = result.crossEntropy(1);
            expect(l2.isFinite());
        }
        /* Overflow test */ {
            const VectorND<float32> result{555.321167, 364.9577942, 355.3863831, -594.8062134};
            for (size_t i = 0; i < result.getLength(); ++i) {
                const float32 s = result.softmax(i);
                expect(s.isFinite());
            }
        }
        /* Diff test */ {
            using dfloat = Diff<float32, DiffMode::Reverse, 1>;
            constexpr int Label = 0;
            const auto x = VectorND<dfloat>::random_uniform<RandomSource>(8);
            auto loss = [&x](size_t k) -> float32 {
                float32 result = 0;
                for (size_t i = 0; i < x.getLength(); ++i)
                    result += x.values().softmax(i) * (float32(i == k) - float32(Label == k));
                return result;
            };

            x.crossEntropy(Label).reverse();
            for (size_t i = 0; i < x.getLength(); ++i)
                expect(scalarNear(x.grads()[i], loss(i), 1E-6));
        }
    }

    void softmaxTest() {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Reverse>;
        const auto factors = VectorND<T>::random_uniform<RandomSource>(8);
        const auto x = VectorND<dfloat>::random_uniform<RandomSource>(8);
        const auto x1 = x;
        VectorND<T> y = softmax(x1.values());
        softmax(x1).reverse(y, factors);

        for (size_t i = 0; i < x.getLength(); ++i)
            x.softmax(i).reverse(factors[i]);
        expect(vectorNear(x.grads(), x1.grads(), 1E-6));
    }

    void testConverts() {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Forward, 1>;
        const VectorND<T> x{0};
        VectorND<dfloat> y = x - T(0.5);
        y[0].grad() = T(0.3);

        VectorND<T> a = y.values();
        expect(a[0] == T(-0.5));

        a = y.grads();
        expect(a[0] == T(0.3));
    }
}

int main() {
    rangeTest<Vector4D<float32>>();
    rangeTest<VectorND<float64>>();
    rangeTest<decltype(std::declval<VectorND<float64>>().tail(0))>();
    strucBindTest();

    crossProductTest();
    innerDotTest();
    hdfTest();
    lnSumExpTest();
    crossEntropyTest();
    softmaxTest();
    testConverts();
    return 0;
}
