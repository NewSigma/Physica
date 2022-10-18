/*
 * Copyright 2022 WeiBo He.
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
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Cycler.h"

using namespace Physica::Core;
using namespace Physica::Core::Parallel;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using ForceModel = Q_TIP4P<ScalarType, PosScalarType>;
constexpr size_t numReplica = 192;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(100);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);

MDCell<ScalarType, PosScalarType> makeSystem(unsigned int cellSize) {
    using LatticeMatrix = typename CrystalCell::LatticeMatrix;
    using PositionMatrix = typename CrystalCell::PositionMatrix;
    const LatticeMatrix lattice{
        -4.6635062604325164,   -0.2499522611778955,    0.0000000000000000,
        -2.1629745970109657,   -4.1943944839773311,    0.0000000000000000,
        -0.2750800827878018,   -0.4169789280520980,   18.0000000000000000
    };
    const PositionMatrix pos {
        0.4553508091084409,  0.3980437584135783,  0.1240303800896787,
        0.4937103263031835,  0.4030549988960055,  0.9488679230950712,
        0.5596918259357793,  0.8517822319914985,  0.1226285591691945,
        0.3686253245184842,  0.6403194088783717,  0.0163388989450929,
        0.5980496296529945,  0.8452761470000689,  0.9474667297556534,
        0.1076134753395919,  0.7454143691003363,  0.1235631952109221,
        0.6847635746074375,  0.9781822048654215,  0.0551553884196571,
        0.9457728898525665,  0.8447981726837335,  0.9479330783198919,
        0.6786065180112657,  0.9738018598142685,  0.1107906720462844,
        0.3747794285285615,  0.6414996922536187,  0.9607035572233092,
        0.7021261874138659,  0.9803177844871507,  0.9573706719037773,
        0.3512600170478342,  0.6392493670479714,  0.1141244566914832
    };

    CrystalCell cell(lattice, pos, {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8}, CrystalCell::Type::Direct);
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    cell.unitToSuper(cellSize, cellSize, 1);
    cell.toCartesian();
    return MDCell<ScalarType, PosScalarType>(std::move(cell));
}

int main() {
    Cycler::init();
    std::mt19937 gen{};

    auto cell = makeSystem(5);
    ForceModel::sortPosition(cell);
    RPMD<ScalarType, PosScalarType> rpmd(std::move(cell), numReplica, 16, temperatureT, thermostatTime, timeStep);
    rpmd.initMomentum(gen);
    ForceModel model(rpmd.phaseToCell(0), pair_cutoff);
    {
        ThreadPool::initThreadPool(4);
        ThreadPool& pool = ThreadPool::getInstance();
        //ProfilerStart("profiler.dat");
        const auto from = Cycler::tic();
        rpmd.nvt_step_for<decltype(gen), decltype(model), ThreadExecutor>(PhyConst<AU>::secondToTime(1 * 1E-13), gen, model);
        const auto to = Cycler::toc();
        //ProfilerStop();
        std::cout << "4 Threads time use: " << Cycler::toSeconds(to - from) << '\n';
        pool.shouldExit();
        ThreadPool::deInitThreadPool();
    }
    return 0;
}
