/*
 * Copyright 2024 Weibo He.
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
#pragma once

#include <fstream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica::Core {
    class MatrixMarket {
        using This = MatrixMarket;
    public:
        ~MatrixMarket() = default;
        /* Static members */
        template<LMatrix T>
        static void read(const char* path, T& target);
    private:
        MatrixMarket();
        MatrixMarket(const This&) = default;
        MatrixMarket(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<LMatrix T>
    void MatrixMarket::read(const char* path, T& target) {
        using ScalarType = typename T::ScalarType;
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
