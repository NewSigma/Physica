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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Utils/Unix/TempFile.h"

using namespace Physica;
using RandomType = Random<MT19937, std::mt19937::default_seed>;

static void crossProductTest() {
    using T = float32;
    VectorND<T> v1{3.845971, 0.000000, 0.000000};
    VectorND<T> v2{-0.007733, 3.835502, 0.000000};
    VectorND<T> v3(v1.crossProduct(v2));
    if (!scalarNear(v3.norm() / T(2), T(7.375614), 1E-7))
        exit(EXIT_FAILURE);
}

static void hdfTest() {
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

static void lnSumExpTest() {
    /* Complex overflow test */ {
        VectorND<cfloat64> v{-1071, -739};
        if (!v.lnSumExp().isFinite())
            exit(EXIT_FAILURE);
    }
    /* Diff test */ {
        using dfloat = Diff<float32, DiffMode::Reverse, 1>;
        const auto x = VectorND<dfloat>::random_uniform<RandomType>(8);
        x.lnSumExp().reverse();

        for (size_t i = 0; i < x.getLength(); ++i) {
            if (!scalarNear(x.grads()[i], x.values().softmax(i), 1E-6))
                exit(EXIT_FAILURE);
        }
    }
}

static void crossEntropyTest() {
    /* Select test */ {
        const VectorND<float32> result{-3.34036088, -109.5531235, 13.51656151, 11.29175949};
        const float32 l1 = result.crossEntropy(3);
        if (!l1.isFinite())
            exit(EXIT_FAILURE);

        const float32 l2 = result.crossEntropy(1);
        if (!l2.isFinite())
            exit(EXIT_FAILURE);
    }
    /* Overflow test */ {
        const VectorND<float32> result{555.321167, 364.9577942, 355.3863831, -594.8062134};
        for (size_t i = 0; i < result.getLength(); ++i) {
            const float32 s = result.softmax(i);
            if (!s.isFinite())
                exit(EXIT_FAILURE);
        }
    }
    /* Diff test */ {
        using dfloat = Diff<float32, DiffMode::Reverse, 1>;
        constexpr int Label = 0;
        const auto x = VectorND<dfloat>::random_uniform<RandomType>(8);
        auto loss = [&x](size_t k) -> float32 {
            float32 result = 0;
            for (size_t i = 0; i < x.getLength(); ++i)
                result += x.values().softmax(i) * (float32(i == k) - float32(Label == k));
            return result;
        };

        x.crossEntropy(Label).reverse();
        for (size_t i = 0; i < x.getLength(); ++i) {
            if (!scalarNear(x.grads()[i], loss(i), 1E-6))
                exit(EXIT_FAILURE);
        }
    }
}

static void softmaxTest() {
    using dfloat = Diff<float32, DiffMode::Reverse>;
    const auto factors = VectorND<float32>::random_uniform<RandomType>(8);
    const auto x = VectorND<dfloat>::random_uniform<RandomType>(8);
    const auto x1 = x;
    VectorND<float32> y = softmax(x1.values());
    softmax(x1).reverse(y, factors);

    for (size_t i = 0; i < x.getLength(); ++i)
        x.softmax(i).reverse(factors[i]);

    if (!vectorNear(x.grads(), x1.grads(), 1E-6))
        exit(EXIT_FAILURE);
}

int main() {
    crossProductTest();
    hdfTest();
    lnSumExpTest();
    crossEntropyTest();
    softmaxTest();
    return 0;
}
