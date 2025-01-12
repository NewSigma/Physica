/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/Schur.h"
#include "Physica/Core/Scalar/Complex.h"

using namespace Physica;

template<Matrix T>
bool isUpperQuasiTriangle(const T& m) {
    if (m.getRow() != m.getCol())
        return false;
    for (size_t i = 0; i < m.getRow() - 1; ++i) {
        if (m(i + 1, i).isZero()) {
            auto col = m.col(i);
            if (!col.tail(i + 1).isZeros())
                return false;
        }
        else if(i < m.getRow() - 2) {
            auto col1 = m.col(i);
            auto col2 = m.col(i + 1);
            if (!col1.tail(i + 2).isZeros() || !col2.tail(i + 2).isZeros())
                return false;
            ++i;
        }
    }
    return true;
}

template<Matrix T>
bool isUpperTriangle(const T& m) {
    if (m.getRow() != m.getCol())
        return false;
    for (size_t i = 0; i < m.getRow() - 1; ++i) {
        auto col = m.col(i);
        if (!col.tail(i + 1).isZeros())
            return false;
    }
    return true;
}

template<Matrix T>
bool realSchurTest(const T& mat, double precision) {
    Schur<typename T::ScalarType> schur(mat, true);
    if (!isUpperQuasiTriangle(schur.getMatrixT()))
        return false;
    T A = schur.getMatrixU() * (schur.getMatrixT() * schur.getMatrixU().transpose()).compute();
    return (A - mat).norm1() <= mat.norm1() * precision;
}

template<Matrix T>
bool schurTest(const T& mat, double precision) {
    static_assert(T::isComplex, "[Error]: Use realSchurTest is prefered");
    Schur<typename T::ScalarType> schur(mat, true);
    if (!isUpperTriangle(schur.getMatrixT()))
        return false;
    T A = schur.getMatrixU() * (schur.getMatrixT() * schur.getMatrixU().hermite()).compute();
    return (A - mat).norm1() <= mat.norm1() * precision;
}

int main() {
    using RealType = float64;
    using ComplexType = Complex<RealType>;
    {
        using MatrixType = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector>;
        const MatrixType mat{{1, 2}, {3, 4}};
        if (!realSchurTest(mat, 1E-15))
            return 1;
    }
    {
        using MatrixType = DenseMatrix<ComplexType, MatrixOption::Col | MatrixOption::Vector>;
        const MatrixType mat{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
        if (!schurTest(mat, 1E-15))
            return 1;
    }
    {
        using MatrixType = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector, 3, 3>;
        const MatrixType mat1{{-149, 537, -27}, {-50, 180, -9}, {-154, 546, -25}};
        if (!realSchurTest(mat1, 1E-14))
            return 1;
        const MatrixType mat2{{-0.590316, -2.19514, -2.37463},
                              {-1.25006, -0.297493, 1.40349},
                              {0.517063, -0.956614, -0.920775}};
        if (!realSchurTest(mat2, 1E-14))
            return 1;
    }
    {
        using MatrixType = DenseMatrix<ComplexType, MatrixOption::Col | MatrixOption::Vector, 3, 3>;
        const MatrixType mat1{{{-149, 37}, {537, -126}, {-27, 0}},
                              {{0, -50}, {0, 180}, {-9, 17}},
                              {{12, -154}, {546, 8}, {-25, 9}}};
        if (!schurTest(mat1, 1E-14))
            return 1;
    }
    {
        using MatrixType = DenseMatrix<RealType, MatrixOption::Row | MatrixOption::Vector>;
        const MatrixType mat{{    -51.01383690006933, 6.924908905363546e-06, 7.198932025559497e-09,-1.272828037030942e-10,-1.443364854227616e-10, 5.138822756983425e-09, 4.861712242670409e-06,-1.277597529382988e-07, 1.017976888750943e-09, 8.941250463097134e-11, 3.242667144581444e-09,  7.65812223572007e-11},
                             { 6.924908905363546e-06,  0.002006525953277094,    0.4127416614600888,   0.02490964847741966,-2.771538634330499e-07,-4.484994220048364e-07,-0.0004996107522897938,  -0.09645908706481159, 0.0007476472528712804, -0.003599700059967799,   -0.1311225590434474, -0.003097747778273297},
                             { 7.198932025559497e-09,    0.4127416614600888,     93.41683282000623,     40.48827018429925,      6.06268744786774, 0.0007298507600543566,    0.6993685886441052,     11.35175309941103,  -0.08855627490432495,    -1.083254145986432,    -39.45448170920827,   -0.9321542370376518},
                             {-1.272828037030942e-10,   0.02490964847741966,     40.48827018429925,     83.11123301592367,     81.67438894071813,  0.002427781007455736,     2.364715877248293,     91.00041086776419,   -0.7063946132535009,    0.6113157992231695,     22.27545753248296,    0.5268876296522678},
                             {-1.443364854227616e-10,-2.771538634330499e-07,      6.06268744786774,     81.67438894071813,    -45.51578644225164,  0.001465204406433001,     1.417894117180667,     55.14028940229678,   -0.4289480169127738,    -2.078228665662546,    -75.69033937713017,    -1.788944253152537},
                             { 5.138822756983425e-09,-4.484994220048364e-07, 0.0007298507600543566,  0.002427781007455736,  0.001465204406433001,    -51.01377984261782,   0.05398341255728714,  0.002424437010296131,-1.686554733501639e-05,-1.120358678570584e-07,-3.407483848815609e-06,-7.663495077598609e-08},
                             { 4.861712242670409e-06,-0.0004996107522897938,    0.6993685886441052,     2.364715877248293,     1.417894117180667,   0.05398341255728714,   0.06123123476312162,     2.402514409700294,  -0.01865575956877191,-0.0002221516784015098, -0.007862265988011478,-0.0001779849155034141},
                             {-1.277597529382988e-07,  -0.09645908706481159,     11.35175309941103,     91.00041086776419,     55.14028940229678,  0.002424437010296131,     2.402514409700294,     97.13038695452103,   -0.7541730794276476,    0.1287911205595803,     4.730168346968156,    0.1119868370907421},
                             { 1.017976888750943e-09, 0.0007476472528712804,  -0.08855627490432495,   -0.7063946132535009,   -0.4289480169127738,-1.686554733501639e-05,  -0.01865575956877191,   -0.7541730794276476,  0.006060106663180818,   -0.1180314801325573,  -0.03351033843348594,  -0.02848921794527029},
                             { 8.941250463097134e-11, -0.003599700059967799,    -1.083254145986432,    0.6113157992231695,    -2.078228665662546,-1.120358678570584e-07,-0.0002221516784015098,    0.1287911205595803,   -0.1180314801325573,      50.2560006099534,    -4.072439521803009,      41.0073828040026},
                             { 3.242667144581444e-09,   -0.1311225590434474,    -39.45448170920827,     22.27545753248296,    -75.69033937713017,-3.407483848815609e-06, -0.007862265988011478,     4.730168346968156,  -0.03351033843348594,    -4.072439521803009,    -62.61826197626559,     -1.79645092586206},
                             {  7.65812223572007e-11, -0.003097747778273297,   -0.9321542370376518,    0.5268876296522678,    -1.788944253152537,-7.663495077598609e-08,-0.0001779849155034141,    0.1119868370907421,  -0.02848921794527029,      41.0073828040026,     -1.79645092586206,    -34.41047976018166}};
        if (!realSchurTest(mat, 1E-10))
            return 1;
    }
    /* Test degeneracy */ {
        using MatrixType = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector, 6, 6>;
        const MatrixType mat1{{ 0.1343184046,             0,             0, -0.1343184056,             0,             0},
                              {            0,  0.1341424528,             0,             0, -0.1341424541,             0},
                              {            0,             0,  0.1342191829,             0,             0, -0.1342191848},
                              {-0.1343184056,             0,             0,  0.1343184065,             0,             0},
                              {            0, -0.1341424541,             0,             0,  0.1341424554,             0},
                              {            0,             0, -0.1342191848,             0,             0,  0.1342191868}};
        if (!realSchurTest(mat1, 1E-15))
            return 1;

        using ComplexMatrix = DenseMatrix<ComplexType, MatrixOption::Col | MatrixOption::Vector, 6, 6>;
        const ComplexMatrix mat2 = mat1;
        if (!schurTest(mat2, 1E-14))
            return 1;
    }
    return 0;
}
