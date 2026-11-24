/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <fstream>
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica {
    class MatrixMarket {
        using This = MatrixMarket;
    public:
        ~MatrixMarket() = default;
        /* Static members */
        static void read(const char* path, Matrix auto& target);
    private:
        MatrixMarket();
        MatrixMarket(const This&) = default;
        MatrixMarket(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    void MatrixMarket::read(const char* path, Matrix auto& target) {
        using ScalarType = std::remove_cvref_t<decltype(target)>::ScalarType;
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: No file found");
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        size_t row, col, numElem;
        fin >> row >> col >> numElem;
        target.resize(row, col);
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        for (size_t i = 0; i < numElem; ++i) {
            size_t r, c;
            ScalarType elem;
            fin >> r >> c >> elem;
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            target(r - 1, c - 1) = elem;
        }

        if (!fin)
            throw BadFileFormatException("[Error]: Bad MTX file");
    }
}
