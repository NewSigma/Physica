/*
 * Copyright 2023-2024 Weibo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <fstream>
#include "Physica/Core/IO/VASP/Poscar.h"
#include "Physica/Core/Physics/MD/ForceModel/QEModel.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/KineticModel/CFireModel.h"
#include "Physica/Core/Physics/MD/Barostat/Berendsen.h"

using namespace Physica::Core;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomType = Random<MT19937, 10002>;
using KineticModel = FreeModel<ScalarType, 3, 1, RPMDIntegrator::Exact>;
using ForceModel = QEModel<ScalarType>;
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;

const static char* input = "&CONTROL\n"
                           "pseudo_dir = '/home/sigma/Program/pot-QE/pslibrary-master/pbe/PSEUDOPOTENTIALS'\n"
                           "calculation = 'scf'\n"
                           "prefix= 'f-BH'\n"
                           "outdir= './output'\n"
                           "tprnfor = .true.\n"
                           "tstress = .true.\n"
                           "/\n"
                           "\n"
                           "&SYSTEM\n"
                           "nat = 12\n"
                           "ntyp = 2\n"
                           "ibrav = 0\n"
                           "ecutwfc = 60\n"
                           "ecutrho = 480\n"
                           "occupations = 'fixed'\n"
                           "/\n"
                           "\n"
                           "&ELECTRONS\n"
                           "mixing_beta = 0.7\n"
                           "conv_thr = 1.0d-10\n"
                           "startingpot = 'file'\n"
                           "startingwfc = 'file'\n"
                           "/\n"
                           "\n"
                           "ATOMIC_SPECIES\n"
                           "H   1.00798  H.pbe-kjpaw_psl.1.0.0.UPF\n"
                           "O   15.9994  O.pbe-nl-kjpaw_psl.1.0.0.UPF\n"
                           "\n"
                           "K_POINTS automatic\n"
                           "4 4 1 0 0 0\n";

MDCell<ScalarType> makeSystem() {
    using CrystalCellType = CrystalCell<ScalarType>;
    using LatticeMatrix = MDCell<ScalarType>::LatticeMatrix;
    using PositionMatrix = MDCell<ScalarType>::PositionMatrix;
    const LatticeMatrix lattice{
         4.6635062604325164,    0.2499522611778955,    0.0000000000000000,
         2.1629745970109657,    4.1943944839773311,    0.0000000000000000,
         0.0000000000000000,    0.0000000000000000,   18.0000000000000000
    };
    PositionMatrix pos {
        0.4553508091084409, 0.3980437584135783, 0.6240303800896787,
        0.4937103263031835, 0.4030549988960055, 0.4488679230950712,
        0.5596918259357793, 0.8517822319914985, 0.6226285591691945,
        0.3686253245184842, 0.6403194088783717, 0.5163388989450929,
        0.5980496296529945, 0.8452761470000689, 0.4474667297556534,
        0.1076134753395919, 0.7454143691003363, 0.6235631952109221,
        0.6847635746074375, 0.9781822048654215, 0.5551553884196571,
        0.9457728898525665, 0.8447981726837335, 0.4479330783198919,
        0.6786065180112657, 0.9738018598142685, 0.6107906720462844,
        0.3747794285285615, 0.6414996922536187, 0.4607035572233092,
        0.7021261874138659, 0.9803177844871507, 0.4573706719037773,
        0.3512600170478342, 0.6392493670479714, 0.6141244566914832
    };
    
    CrystalCellType cell({lattice, pos, CrystalCellType::Type::Direct}, {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8});
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    MDCell<ScalarType> cell1(std::move(cell));
    //for (size_t i = 1; i < cell1.getNumParticle(); ++i)
    //    cell1.setMass(i, cell1.getMass(0)); // Set all mass equal to mass H is vital to performance
    pos = cell1.getPos();
    std::normal_distribution<double> dist(0, 0.5);
    for (auto& elem : pos.asArray())
        elem += ScalarType::random_any<decltype(dist), RandomType>(dist); //Perturbation
    cell1.setPos(pos);
    return cell1;
}

ForceModel makeForceModel() {
    Poscar<ScalarType>::ElementTypeArray elementArray{1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8};
    TempFile inputFile("/tmp/tmpXXXXXX");
    std::ofstream os(inputFile.getName());
    os << input;
    os.close();
    return ForceModel("/home/sigma/Program/qe-7.0/bin/pw.x", inputFile.getName(), std::move(elementArray), 4);
}

int main() {
    RPMD<ScalarType, 3, 1> rpmd(makeSystem(), 1, 1, 0, timeStep);

    ForceModel forceModel = makeForceModel();
    KineticModel kineticModel(0, 1);
    FireModel<ScalarType, 3> fire(timeStep, 10 * timeStep);

    VectorType f2norm{};
    f2norm.reserve(300);
    for (size_t i = 0; i < 300; ++i) {
        rpmd.fire_vstep<KineticModel, ForceModel, SequentialExecutor>(fire, kineticModel, forceModel);
        f2norm.append(fire.getForceNorm());
        const ScalarType maxForce = abs_elem(rpmd.getForce()).max();
        if (maxForce < ScalarType(1E-5 * PhyConst<QE>::planck / PhyConst<AU>::planck))
            break;
    }
    H5File h5f("force.h5", H5File::ReadWrite | H5File::Creat);
    f2norm.write(h5f, "force");
    rpmd.phaseToCell(0).write(h5f, "cell");
    return 0;
}
