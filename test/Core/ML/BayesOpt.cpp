/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/ML/BayesOpt.h"
#include "Physica/Core/ML/Optimizer/Adadelta.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Optimization/TestFunction.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using dfloat = Diff<T, DiffMode::Reverse>;
using RandomSource = Random<>;

namespace {
    void basic() {
        auto fn = [](const VectorND<T>& x) static noexcept {
            return -forrester(x.front());
        };

        GPModel<dfloat> model(1);
        auto x = VectorND<T>::linspace(0, 1, 4);
        auto y = VectorND<T>::generate([&](size_t i) noexcept {
            return fn({x[i]});
        }, x.getLength());

        using Cube = BayesOpt<T>::Cube;
        const size_t argmax = y.argmax();
        auto bayes = BayesOpt<T>(Vegas<T, true>(Cube{0, 1}, 1000, 30), {x[argmax]}, y[argmax]);
        auto opt = Adadelta<T>{};
        for (int i = 0; i < 10; ++i) {
            for (int _ = 0; _ < 128; ++_) {
                model.regression(x.transpose(), y, T(std::numeric_limits<T>::epsilon())).reverse(-1);
                model.step(opt);
                model.zero_grad();
            }

            auto [xnew, ynew] = bayes.template propose<RandomSource>(fn, model, 64);
            x.append(xnew.front());
            y.append(ynew);
        }
        expect(fn(bayes.getOptimal()) == bayes.getMaximum());
        expect<RandomSource>(scalarNear(bayes.getOptimal().front(), T(0.75724875784183199), 1E-2));
    }
}

int main() {
    basic();
    return 0;
}
