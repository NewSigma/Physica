/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FDM_H
#define PHYSICA_FDM_H

#include <QtCore/QVector>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/ColumnMatrix.h"
/*!
 * Solve the laplace equation using FDM. This is created for study and is very shabby.
 */
namespace Physica::Core {
    class FDM {
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
        static const int iterateMax;
        static const double precision;
        ColumnMatrix data;
        QVector<Boundary> boundaries;
    public:
        FDM(size_t column, size_t row);

        void addBoundary(const Boundary& boundary, const Scalar<MultiPrecision, false>& value);
        void loop();
        /* Getters */
        [[nodiscard]] const ColumnMatrix& getData() { return data; }
    private:
        bool onBoundary(size_t column, size_t row);
    };
}

#endif