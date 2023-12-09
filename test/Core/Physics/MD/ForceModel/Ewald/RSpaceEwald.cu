/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Core/Math/Random/RandomPool.h"

using namespace Physica::Core;
using ScalarType = Scalar<Float>;
using MDCellType = MDCell<ScalarType>;
using LatticeMatrix = typename MDCellType::LatticeMatrix;
using PositionMatrix = typename MDCellType::PositionMatrix;
using MassVector = typename MDCellType::MassVector;
using HostForceModel = RSpaceEwald<ScalarType>;
using DeviceForceModel = device_obj<HostForceModel>;
using RandomPoolType = RandomPool<std::mt19937, 10000>;

template<class RandomGenerator>
MDCellType makeSystem(size_t numMolecular, RandomGenerator& gen) {
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    pos.random_uniform(gen);
    typename MDCellType::MassVector massVec(numMolecular, 1.0);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));
    cell.scale(ScalarType(20));
    return cell;
}

int main() {
    auto& gen = RandomPoolType::getGen();
    const auto cell = makeSystem(108, gen);
    const auto& pos = cell.getPos();
    Vector<ScalarType> charges(cell.getNumParticle(), 1.0);
    auto tail = charges.tail(cell.getNumParticle() / 2);
    tail = ScalarType(-1);

    HostForceModel hostModel(cell.getLattice(), charges);
    DeviceForceModel deviceModel(cell.getLattice(), charges);

    const auto f0 = hostModel.template force_short<SequentialExecutor>(pos);
    const auto f1 = deviceModel.template force_short<CudaExecutor>(pos);
    if (!vectorNear(f0, f1, 1E-4))
        return 1;
    return 0;
}
