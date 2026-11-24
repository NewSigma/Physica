/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/ManyBody/GreenSampler/ImagMagnon.h"
#include "Test.h"

using namespace Physica;

int main() {
    syntax_only([]() {
        using T = cfloat64;
        ImagMagnon<T> sampler({}, 1UZ);
        std::ignore = sampler.getNumSite();
        std::ignore = sampler.getNumSplit();
        std::ignore = sampler.calcMean();
        std::ignore = sampler.calcSign();
    });
    return 0;
}
