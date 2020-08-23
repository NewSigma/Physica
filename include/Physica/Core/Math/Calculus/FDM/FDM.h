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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
/*!
 * Solve the laplace equation using FDM. This is created for study and is very shabby.
 */
namespace Physica::Core {
    class FDMBase {
    public:
        static const int iterateMax;
        static const double precision;
    };

    template<class T>
    class FDM : public FDMBase {
    public:
        enum BoundaryType {
            Row,
            Column
        };

        struct Boundary {
            BoundaryType type;
            size_t index, from, to;
        };
    private:
        Matrix<T, Column, Dynamic, Dynamic> data;
        QVector<Boundary> boundaries;
    public:
        FDM(size_t column, size_t row);

        void addBoundary(const Boundary& boundary, const T& value);
        void loop();
        /* Getters */
        [[nodiscard]] const Matrix<T, Column, Dynamic, Dynamic>& getData() { return data; }
    private:
        bool onBoundary(size_t column, size_t row);
    };
}

#include "FDMImpl.h"

#endif