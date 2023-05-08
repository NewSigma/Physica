/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Flatten.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/Ewald.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

using namespace Physica;
using namespace Physica::Core;
using namespace Physica::Utils;

template<class ScalarType>
void VASPTest() {
    constexpr static bool isFloat = ScalarType::option == Float;
    constexpr double prec = isFloat ? 2E-5 : 1E-5;

    const double lengthInBohr = PhyConst<AU>::angstormToBohr(3);
    CrystalCell cell({{lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {0.5, 0.5, 0.5}, CrystalCell::Type::Direct}, {14});
    Ewald<ScalarType> ewald(cell.getLattice(), {4});
    const auto energy = ewald.potentialEnergy(cell.getPos());
    if (!scalarNear(energy, ScalarType(PhyConst<AU>::eVToHartree(-108.95061336198556)), prec))
        exit(EXIT_FAILURE);
}
/**
 * Reference:
 * [1] pyewald(github.com/lukeolson/pyewald)
 */
template<class ScalarType>
void madelungTest() {
    constexpr static bool isFloat = ScalarType::option == Float;
    {
        const double lengthInBohr = PhyConst<AU>::angstormToBohr(5.6903014761756712);
        CrystalCell NaCl({{lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {
                            0.0, 0.0, 0.0,
                            0.0, 0.5, 0.5,
                            0.5, 0.0, 0.5,
                            0.5, 0.5, 0.0,
                            0.5, 0.5, 0.5,
                            0.5, 0.0, 0.0,
                            0.0, 0.5, 0.0,
                            0.0, 0.0, 0.5
                        }, CrystalCell::Type::Direct}, {1, 1, 1, 1, 1, 1, 1, 1});
        NaCl.toCartesian();
        Ewald<ScalarType> ewald(NaCl.getLattice(), {1, 1, 1, 1, -1, -1, -1, -1});
        const auto energy = ewald.potentialEnergy(NaCl.getPos());
        const auto madelung = -(energy / 4) * (lengthInBohr / 2); //We have 4x unit cell so energy is divided by 4
        constexpr double prec = isFloat ? 1E-6 : 1E-7;
        if (!scalarNear(madelung, ScalarType(1.7475645946331822), prec))
            exit(EXIT_FAILURE);
    }
    {
        const double lengthInBohr = 1;
        CrystalCell CsCl({{lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {
                            0.0, 0.0, 0.0,
                            0.5, 0.5, 0.5,
                         }, CrystalCell::Type::Cartesian}, {1, 1});
        Ewald<ScalarType> ewald(CsCl.getLattice(), {1, -1});
        const auto energy = ewald.potentialEnergy(CsCl.getPos());
        const auto madelung = -energy * (lengthInBohr * 0.5 * std::sqrt(3.0));
        constexpr double prec = isFloat ? 1E-6 : 1E-9;
        if (!scalarNear(madelung, ScalarType(1.76267477307099), prec))
            exit(EXIT_FAILURE);
    }
    {
        const double lengthInBohr = 0.5;
        CrystalCell ZnS({{0, lengthInBohr, lengthInBohr, lengthInBohr, 0, lengthInBohr, lengthInBohr, lengthInBohr, 0}, {
                            0.0, 0.0, 0.0,
                            0.25, 0.25, 0.25,
                        }, CrystalCell::Type::Direct}, {1, 1});
        ZnS.toCartesian();
        Ewald<ScalarType> ewald(ZnS.getLattice(), {1, -1});
        const auto energy = ewald.potentialEnergy(ZnS.getPos());
        const auto madelung = -energy * (lengthInBohr * 0.5 * std::sqrt(3.0));
        constexpr double prec = isFloat ? 1E-6 : 1E-9;
        if (!scalarNear(madelung, ScalarType(1.63805505338879), prec))
            exit(EXIT_FAILURE);
    }
}

template<class ScalarType>
void forceTest() {
    constexpr static bool isFloat = ScalarType::option == Float;
    constexpr double prec = isFloat ? 2E-1 : 1E-3; //Poor precision at single precision

    using Executor = SequentialExecutor;
    using LatticeMatrix = typename CrystalCell::LatticeMatrix;
    using PositionMatrix = typename CrystalCell::PositionMatrix;
    using ScalarType_ = typename CrystalCell::ScalarType;
    LatticeMatrix lattice{
        4.6635062604325164,   0.2499522611778955,    0.0000000000000000,
        2.1629745970109657,   4.1943944839773311,    0.0000000000000000,
        0.2750800827878018,   0.4169789280520980,   18.0000000000000000
    };
    PositionMatrix pos {
        3.018608093,  1.835086465,  2.232546806, 
        3.435233593,  2.209633827,  17.07962227, 
        4.486242294,  3.763740540,  2.207314014,
        3.108575344,  2.784703970,  0.294100195,
        4.877948284,  4.089979172,  17.05440140,
        2.148158073,  3.204983473,  2.224137545, 
        5.324354172,  4.297039032,  0.992797017,
        6.498651505,  4.175083160,  17.06279564,
        5.301470280,  4.300325871,  1.994232059,
        3.399604082,  3.184972763,  17.29266357,
        5.658125877,  4.686541080,  17.23267174,
        3.052176714,  2.816649675,  2.054240227
    };
    lattice *= reciprocal(ScalarType_(PhyConst<SI>::bohrRadius * 1E10));
    pos *= reciprocal(ScalarType_(PhyConst<SI>::bohrRadius * 1E10));
    const Ewald<ScalarType, ScalarType_> ewald(lattice, {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8});
    const auto force1 = ewald.template force<Executor>(pos);
    Vector<ScalarType> force2(force1.getLength());
    /* Force from finite differential */ {
        for (size_t i = 0; i < force2.getLength(); ++i) {
            force2[i] = -Differential<ScalarType>::ridders([i, &pos, &ewald](ScalarType x) -> ScalarType {
                PositionMatrix temp = pos;
                *(temp.begin() + i) = ScalarType_(x);
                return ewald.potentialEnergy(temp);
            }, pos.flatten().calc(i), 0.3);
        }
    }

    if (!vectorNear(force1, force2, prec))
        exit(EXIT_FAILURE);
}

int main() {
    VASPTest<Scalar<Double, false>>();
    VASPTest<Scalar<Float, false>>();
    madelungTest<Scalar<Double, false>>();
    madelungTest<Scalar<Float, false>>();
    forceTest<Scalar<Double, false>>();
    forceTest<Scalar<Float, false>>();
    return 0;
}
