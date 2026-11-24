/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2019-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Real.h"
#include "Test.h"

using namespace Physica;

int main() {
    /* test Pi */ {
        auto Pi = MathConst<Real<FloatMP>>::getInstance().PI;
        expect(Pi.getSize() == 5
            && Pi[0] == 11424456171093639400UL
            && Pi[1] == 11820040416388919749UL
            && Pi[2] == 1376283091369227076
            && Pi[3] == 2611923443488327891
            && Pi[4] == 3);
    }
    /* test E */ {
        auto E = MathConst<Real<FloatMP>>::getInstance().E;
        expect(E.getSize() == 4
            && E[0] == 7126689189968796226
            && E[1] == 13794904443024896967UL
            && E[2] == 13249961062380153450UL
            && E[3] == 2);
    }
    return 0;
}
