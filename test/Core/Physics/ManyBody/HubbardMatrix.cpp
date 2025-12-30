/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/KFermiRepr.h"
#include "Test.h"

using namespace Physica;
constexpr double HoppingT = 1.0;
constexpr double RepelU = 2;

namespace {
    template<unsigned int NumSite, unsigned int NumSpinUp, unsigned int NumSpinDown>
    void testRSpinMatrix1D() {
        using T = float64;
        using ReprType = FermiRepr<1, NumSite, NumSpinUp == NumSpinDown>;

        SquareLattice<1> lattice({NumSite}, 1);
        const HubbardMatrix<T, ReprType> hamilton(HoppingT, RepelU, lattice, ReprType(NumSpinUp, NumSpinDown));
        const size_t numState = hamilton.getNumState();
        MatrixND<T> mat(numState, numState);
        for (size_t i = 0; i < numState; ++i) {
            auto col = mat.col(i);
            col = hamilton * UnitVector<T>(i, numState);
        }
        expect(matrixNear(hamilton, mat, 1E-15));
    }

    template<unsigned int NumSiteX, unsigned int NumSiteY, unsigned int NumSpinUp, unsigned int NumSpinDown>
    void testRSpinMatrix2D() {
        using T = float64;
        using ReprType = FermiRepr<2, NumSiteX * NumSiteY, NumSpinUp == NumSpinDown>;

        SquareLattice<2> lattice({NumSiteX, NumSiteY}, 1);
        const HubbardMatrix<T, ReprType> hamilton(HoppingT, RepelU, lattice, ReprType(NumSpinUp, NumSpinDown));
        const size_t numState = hamilton.getNumState();
        MatrixND<T> mat(numState, numState);
        for (size_t i = 0; i < numState; ++i) {
            auto col = mat.col(i);
            col = hamilton * UnitVector<T>(i, numState);
        }
        expect(matrixNear(hamilton, mat, 1E-15));
    }

    void testKSpinMatrix() {
        constexpr unsigned int NumSite = 4;
        constexpr unsigned int NumParticle = NumSite / 2;
        using T = cfloat64;
        using ReprType = KFermiRepr<1, NumSite, true>;

        SquareLattice<1> lattice({NumSite}, 1);
        const HubbardMatrix<T, ReprType> hamilton(HoppingT, RepelU, lattice, ReprType({NumParticle, NumParticle}, 0));
        const size_t numState = hamilton.getNumState();
        MatrixND<T> mat(numState, numState);
        for (size_t i = 0; i < numState; ++i) {
            auto col = mat.col(i);
            col = hamilton * UnitVector<T>(i, numState);
        }
        expect(matrixNear(hamilton, mat, 1E-15));
    }

    template<BoundaryCond BC>
    void testVecProduct() {
        constexpr int Dim = 1;
        constexpr int NumSite = 4;
        using T = std::conditional<BC == BoundaryCond::TBC, cfloat64, float64>::type;
        using ReprType = FermiRepr<1, NumSite, false>;
        using Hamilton = HubbardMatrix<T, ReprType, BC>;
        using RandomSource = Random<>;

        SquareLattice<Dim, BC> lattice;
        if constexpr (BC == Physica::BoundaryCond::TBC)
            lattice = SquareLattice<Dim, BC>({NumSite}, 1, {float64::template random_uniform<RandomSource>() * MathConst<float64>::pi});
        else
            lattice = SquareLattice<Dim, BC>({NumSite}, 1);
        const Hamilton hamiltonH(HoppingT, RepelU, lattice, ReprType(2, 1));

        const auto v = VectorND<float64>::random_uniform<RandomSource>(hamiltonH.getNumState());
        const VectorND<T> v1 = hamiltonH * v;
        VectorND<T> v2(v.getLength());
        for (size_t i = 0; i < v.getLength(); ++i)
            v2[i] = (hamiltonH * v).calc(i);
        expect(vectorNear(v1, v2, 1E-14));
    }
}

int main() {
    testRSpinMatrix1D<3, 2, 1>();
    testRSpinMatrix1D<2, 1, 1>();
    testRSpinMatrix2D<2, 3, 2, 3>();
    testKSpinMatrix();
    testVecProduct<BoundaryCond::PBC>();
    testVecProduct<BoundaryCond::TBC>();
    return 0;
}
