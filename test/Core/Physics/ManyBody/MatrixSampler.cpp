/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/ManyBody/GreenSampler/MatrixSampler.h"
#include "Test.h"

using namespace Physica;

namespace {
    void complex_instantiation() {
        syntax_only([]() {
            using T = cfloat64;
            MatrixSampler<T> sampler({}, {}, {}, {});
            sampler.sample({}, {});
            std::ignore = sampler.calcRawMean();
            std::ignore = sampler.calcMean();
            std::ignore = sampler.calcStructFactor();
            std::ignore = sampler.getObserves();
            std::ignore = sampler.getNumSiteX();
            std::ignore = sampler.getNumSiteY();
            std::ignore = sampler.getNumSite();
            std::ignore = sampler.getNumSample();
        });
    }
}

int main() {
    complex_instantiation();
    return 0;
}
