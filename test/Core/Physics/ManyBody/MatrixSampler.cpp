/*
 * Copyright 2025-2026 Weibo He.
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
