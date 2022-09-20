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
#include <iostream>
#include <fstream>
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using namespace Physica::Core::Parallel;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using ForceModel = Q_TIP4P<ScalarType, PosScalarType>;
constexpr size_t numReplica = 32;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.1;
constexpr unsigned int numMolecular = 32;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
constexpr double massMoleculeInSI = PhyConst<SI>::atomMass(1) * 2 + PhyConst<SI>::atomMass(8);
constexpr double cellVolume = ((numMolecular * massMoleculeInSI * 1000 / 0.997) * 1E-6) / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius);

void testHytrogenList() {
    typename CrystalCell::LatticeMatrix lattice{
        4.6635062604325164,   0.2499522611778955,    0.0000000000000000,
        2.1629745970109657,   4.1943944839773311,    0.0000000000000000,
        0.2750800827878018,   0.4169789280520980,   18.0000000000000000
    };
    typename CrystalCell::PositionMatrix pos {
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
    CrystalCell unit{std::move(lattice), std::move(pos), {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8}, CrystalCell::Type::Direct};
    const MDCell<ScalarType, PosScalarType> cell(std::move(unit));
    Q_TIP4P<ScalarType, PosScalarType> model(cell, 1, 1E-4);
    model.guessHytrogenList(cell);
    if (!model.checkHytrogenList())
        exit(EXIT_FAILURE);
}

template<class RandomGenerator>
Vector<PosScalarType, 3> randomVector(RandomGenerator& gen) {
    std::uniform_real_distribution dist{};
    const PosScalarType theta(dist(gen) * M_PI);
    const PosScalarType phi(dist(gen) * M_PI * 2);
    Vector<PosScalarType, 3> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    result *= PosScalarType(ForceModel::equalR);
    return result;
}

template<class RandomGenerator>
RPMD<ScalarType, PosScalarType> makeSystem(RandomGenerator& gen) {
    const size_t maxIndexH = numMolecular * 2;
    const size_t maxIndexO = numMolecular * 3;
    const size_t numAtom = maxIndexO;

    const Scalar<Float, false> latticeFactor(std::cbrt(cellVolume));
    typename CrystalCell::LatticeMatrix lattice = CrystalCell::LatticeMatrix::unitMatrix(3);
    lattice *= latticeFactor;

    typename CrystalCell::PositionMatrix pos(numAtom, 3);
    std::uniform_real_distribution dist{};
    for (size_t i = 0; i < numMolecular; ++i) {
        auto posO = pos.row(i + maxIndexH);
        posO[0] = latticeFactor * dist(gen);
        posO[1] = latticeFactor * dist(gen);
        posO[2] = latticeFactor * dist(gen);
        auto posH1 = pos.row(2 * i);
        auto posH2 = pos.row(2 * i + 1);
        posH1 = posO + randomVector(gen);
        posH2 = posO + randomVector(gen);
    }

    typename CrystalCell::AtomicArray atomicNumbers(numAtom);
    for (size_t i = 0; i < maxIndexH; ++i)
        atomicNumbers[i] = 1;
    for (size_t i = maxIndexH; i < maxIndexO; ++i)
        atomicNumbers[i] = 8;

    CrystalCell cell(std::move(lattice), std::move(pos), std::move(atomicNumbers), CrystalCell::Type::Cartesian);
    return RPMD<ScalarType, PosScalarType>(MDCell<ScalarType, PosScalarType>(std::move(cell)), numReplica, temperatureT, thermostatTime, timeStep);
}

void testMD() {
    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    auto rpmd = makeSystem(gen);
    rpmd.initMomentum(gen);
    ForceModel model(rpmd.phaseToCell(0), pair_cutoff, 1E-4);

    ThreadPool::initThreadPool(4);
    ThreadPool& pool = ThreadPool::getInstance();
    {
        rpmd.updateForce<decltype(model), ThreadExecutor>(model);
        rpmd.nvt_step_for<decltype(gen), decltype(model), ThreadExecutor>(PhyConst<AU>::secondToTime(2 * 1E-12), gen, model);
        ScalarType bond = 0;
        for (size_t i = 0; i < 100; ++i) {
            auto pos = rpmd.getPos();
            ScalarType temp = 0;
            for (size_t j = 0; j < numMolecular; ++i) {
                auto posO = pos.row(2 * numMolecular + j);
                auto posH1 = pos.row(model.getHytrogenList()[j].first);
                auto posH2 = pos.row(model.getHytrogenList()[j].second);
                toNextMean(temp, 2 * j, (posO.asVector() - posH1).norm());
                toNextMean(temp, 2 * j + 1, (posO.asVector() - posH2).norm());
            }
            rpmd.nvt_step<decltype(gen), decltype(model), ThreadExecutor>(gen, model);
            toNextMean(bond, i, temp);
        }
        std::cout << PhyConst<AU>::bohrToAngstorm(double(bond)) << std::endl;
        std::ofstream fout("phase");
        fout << rpmd.getPhasePos();
    }
    pool.shouldExit();
    ThreadPool::deInitThreadPool();
}

int main() {
    testHytrogenList();
    testMD();
    return 0;
}
