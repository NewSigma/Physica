/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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