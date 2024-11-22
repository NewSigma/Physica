/*
 * Copyright 2023-2024 Weibo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <iostream>
#include <fstream>
#include <QApplication>
#include "Physica/Core/IO/VASP/Poscar.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/BKSModel.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/KineticModel/FireModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Utils/Unix/UnixHelper.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomType = Random<MT19937>;
using MDType = RPMD<ScalarType, 3, 1>;
using ThermostatType = Langevin<ScalarType, 3, 1>;
using KineticModel = FreeModel<ScalarType, 3, 1, RPMDIntegrator::Exact>;
using ForceModel = BKSModel<ScalarType, Ewald<ScalarType, RSpaceEwald<ScalarType, true>>, false>;
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(10);
constexpr size_t numStep = 1000;
/**
 * Initial structure from [1], modifying lattice according to [2]
 * Reference:
 * [1] Materials Project, mp-7000; https://doi.org/10.17188/1272685
 * [2] Comput. Mater. Sci. 175, 109584 (2020); https://doi.org/10.1016/j.commatsci.2020.109584
 */
MDCell<ScalarType> makeSystem() {
    using CrystalCellType = CrystalCell<ScalarType>;
    using LatticeMatrix = MDCell<ScalarType>::LatticeMatrix;
    using PositionMatrix = MDCell<ScalarType>::PositionMatrix;
    const LatticeMatrix lattice{
        5.0000000000, 0.0000000000, 0.0000000000,
       -2.5000000000, 4.3333333333, 0.0000000000,
        0.0000000000, 0.0000000000, 5.3333333333
    };
    PositionMatrix pos {
        0.256094009, 0.414853990, 0.794543028,
        0.585146010, 0.841239989, 0.127876326,
        0.158759996, 0.743906021, 0.461209685,
        0.743906021, 0.158759996, 0.538790345,
        0.414853990, 0.256094009, 0.205457002,
        0.841239989, 0.585146010, 0.872123659,
        0.523694992, 0.523694992, 0.000000000,
        0.000000000, 0.476305008, 0.666666687,
        0.476305008, 0.000000000, 0.333333343
    };
    
    CrystalCellType cell({lattice, pos, CrystalCellType::Type::Direct}, {8, 8, 8, 8, 8, 8, 14, 14, 14});
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    MDCell<ScalarType> cell1(std::move(cell));
    cell1.toSuperCell<ExtendCellOption::AtomMajor>(4, 6, 3);
    /* To orthogonal */ {
        LatticeMatrix latt = cell1.getLattice();
        latt(1, 0) = ScalarType(0);
        cell1.setLattice(std::move(latt));
        cell1.normalize();
    }
    return cell1;
}

int main() {
    MDType rpmd(makeSystem(), 1, 1, 0, timeStep);
    const char* oldDataPath = "/kaggle/input/tempsio2/SiO2.h5";
    if (fileExists(oldDataPath)) {
        H5File h5f(oldDataPath, H5File::ReadOnly);
        rpmd.read(h5f, "md");
    }

    VectorType energy(numStep);
    ForceModel forceModel(rpmd.phaseToCell(0), pair_cutoff, {});
    for (size_t step = 0; step <= numStep; ++step) {
        ScalarType temperatureT = ScalarType(PhyConst<AU>::kToTemperature(step * 10));
        
        KineticModel kineticModel(temperatureT, 1);
        ThermostatType thermo(temperatureT, thermostatTime, true);

        rpmd.nvt_step_for<ThermostatType, RandomType, KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(1E-12), thermo, RandomType::getInstance(), kineticModel, forceModel);
        energy[step] = rpmd.calcPotential<ForceModel, SequentialExecutor>(forceModel);
    }

    H5File h5f("SiO2.h5", H5File::ReadWrite | H5File::Creat);
    rpmd.write(h5f, "md");
    energy.write(h5f, "E");

    const size_t numParticle = rpmd.getNumParticle();
    Poscar<ScalarType> poscar(rpmd.makeAverageCell(), {8, 14}, {numParticle * 2 / 3, numParticle / 3});
    std::ofstream fout("POSCAR");
    fout << poscar;
    return 0;
}
