/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_FDM_H
#define PHYSICA_FDM_H

#include <QtCore/QVector>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
/*!
 * Solve the laplace equation using FDM. This is created for study and is very shabby.
 */
namespace Physica::Core {
    class FDMBase {
    public:
        enum BoundaryType {
            Row,
            Column
        };
        /*!
         * Only supports horizontal or vertical Dirichlet boundary conditions.
         */
        struct Boundary {
            BoundaryType type;
            //index of row or column
            size_t index;
            //index of node we start from
            size_t from;
            //index of node we end to
            size_t to;
        };
    public:
        static const int iterateMax;
        static const double precision;
    };
    /*!
     * Solve 2D laplacian function using FDM.
     */
    template<class T>
    class FDM : public FDMBase {
    private:
        DenseMatrix<T, DenseMatrixType::VectorColumn, Dynamic, Dynamic> data;
        QVector<Boundary> boundaries;
    public:
        FDM(size_t column, size_t row);
        /* Operations */
        void addBoundary(const Boundary& boundary, const T& value);
        void loop();
        /* Getters */
        [[nodiscard]] const DenseMatrix<T, DenseMatrixType::VectorColumn, Dynamic, Dynamic>& getData() { return data; }
    private:
        bool onBoundary(size_t column, size_t row);
    };
}

#include "FDMImpl.h"

#endif