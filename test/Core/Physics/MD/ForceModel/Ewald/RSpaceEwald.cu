/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/Ewald/RSpaceEwald.cuh"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using ScalarType = float32;
using MDCellType = MDCell<ScalarType>;
using LatticeMatrix = MDCellType::LatticeMatrix;
using PositionMatrix = MDCellType::PositionMatrix;
using MassVector = MDCellType::MassVector;
using HostForceModel = RSpaceEwald<ScalarType>;
using DeviceForceModel = device_obj<HostForceModel>;
using RandomSource = Random<MT19937, 10000>;

namespace {
    MDCellType makeSystem(size_t numMolecular) {
        MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
        auto pos = MDCellType::PositionMatrix::template random_uniform<RandomSource>(numMolecular, 3);
        MDCellType::MassVector massVec(numMolecular, 1.0);
        MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));
        cell.scale(ScalarType(20));
        return cell;
    }
}

int main() {
    const auto cell = makeSystem(108);
    const auto& pos = cell.getPos();
    VectorND<ScalarType> charges(cell.getNumParticle(), 1.0);
    auto tail = charges.tail(cell.getNumParticle() / 2);
    tail = ScalarType(-1);

    HostForceModel hostModel(cell.getLattice(), charges);
    DeviceForceModel deviceModel(cell.getLattice(), charges);
    {
        const auto f0 = hostModel.template force_short<Sequential>(pos);
        const auto f1 = deviceModel.template force_short<GPU>(pos);
        if (!vectorNear(f0, f1, 1E-4))
            return 1;
    }
    {
        const auto v0 = hostModel.virial(pos);
        const auto v1 = deviceModel.virial(pos);
        if (!matrixNear(v0, v1, 1E-4))
            return 1;
    }
    return 0;
}
