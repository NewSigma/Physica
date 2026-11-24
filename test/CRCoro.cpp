/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/CRCoro.h"
#include "Test.h"

using namespace Physica;

namespace {
    struct A : public CRCoro<A> {
        int i;
        A() = default;
        A(int i_) : i(i_) {}
    };

    A fn() {
        const A x{123};
        co_return x; // Regression test for returning a const value
    }
}

int main() {
    expect(fn().i == 123);
}
