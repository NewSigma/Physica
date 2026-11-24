/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;

namespace {
    void simplify() {
        syntax_only([]() {
            using T = float32;
            auto x = -(VectorND<T>() * T(2));
            static_assert(x.getExprID() == ExprID::Mul);
        });
    }

    void regression() {
        syntax_only([]() {
            std::ignore = VectorND<cfloat32>{} - VectorND<float32>(); // Test that we do not diagnostic this
        });
    }
}

int main() {
    simplify();
    regression();
    return 0;
}
