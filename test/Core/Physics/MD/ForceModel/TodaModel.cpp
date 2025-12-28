/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/MD/ForceModel/TodaModel.h"

using namespace Physica;

using ScalarType = float64;
using RandomSource = Random<MT19937, 10000>;
using MDCellType = MDCell<ScalarType, 1>;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double unitMassM = 1;

namespace {
    MDCellType makeSystem() {
        MDCellType::LatticeMatrix lattice{latticeSize};

        auto posVec = VectorND<ScalarType>::random_uniform<RandomSource>(numMolecular);
        posVec *= ScalarType(latticeSize);
        std::sort(posVec.begin(), posVec.end());
        MDCellType::PositionMatrix pos(numMolecular, 1);
        pos.col(0) = posVec;

        MDCellType::MassVector massVec(numMolecular);
        for (size_t i = 0; i < numMolecular; ++i) {
            massVec[i] = (i % 2U == 0) ? unitMassM : (unitMassM * 10);
        }
        return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
    }

    template<bool IsPeriodBoundary>
    void forceConstTest() {
        using ForceModel = TodaModel<ScalarType, IsPeriodBoundary>;
        ForceModel model(1.0);
        const auto cell = makeSystem();
        const auto fc = model.forceConst(cell);
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < fc.getRow(); ++i)
                if (!scalarNear(fc.row(i).sum(), ScalarType(0), 1E-15))
                    exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < fc.getRow(); ++i) {
            for (size_t j = 0; j < fc.getCol(); ++j) {
                if (!scalarNear(fc[i, j], model.forceConst(cell, i, j), 1E-15))
                    exit(EXIT_FAILURE);
            }
        }
    }
}

int main() {
    forceConstTest<false>();
    forceConstTest<true>();
    return 0;
}
