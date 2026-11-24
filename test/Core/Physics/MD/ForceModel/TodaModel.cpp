/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/MD/ForceModel/TodaModel.h"
#include "Test.h"

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
                expect(scalarNear(fc.row(i).sum(), ScalarType(0), 1E-15));
        }

        for (size_t i = 0; i < fc.getRow(); ++i)
            for (size_t j = 0; j < fc.getCol(); ++j)
                expect(scalarNear(fc[i, j], model.forceConst(cell, i, j), 1E-15));
    }
}

int main() {
    forceConstTest<false>();
    forceConstTest<true>();
    return 0;
}
