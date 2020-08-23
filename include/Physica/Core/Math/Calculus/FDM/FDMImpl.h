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
#ifndef PHYSICA_FDMIMPL_H
#define PHYSICA_FDMIMPL_H

namespace Physica::Core {
    template<class T>
    FDM<T>::FDM(size_t column, size_t row) : data(column) {
        for(size_t i = 0; i < row; ++i)
            data.allocate(Vector<MultiScalar>::zeroVector(row), i);
    }
    /*!
     * By default, edge of the matrix is set zero.
     */
    template<class T>
    void FDM<T>::addBoundary(const Boundary& boundary, const T& value) {
        boundaries.push_back(boundary);
        if(boundary.type == Row) {
            for(size_t i = boundary.from; i <= boundary.to; ++i)
                data(boundary.index, i) = value;
        }
        else {
            for(size_t i = boundary.from; i <= boundary.to; ++i)
                data(i, boundary.index) = value;
        }
    }

    template<class T>
    void FDM<T>::loop() {
        const auto column_1 = data.column() - 1;
        const auto row_1 = data.row() - 1;
        int iterate = 0;
        MultiScalar copy;
        const MultiScalar limit(precision);
        bool keep;
        do {
            keep = false;
            for(size_t i = 1; i < column_1; ++i) {
                for(size_t j = 1; j < row_1; ++j) {
                    if(!onBoundary(i, j)) {
                        auto& temp = data(j, i);
                        copy = std::move(temp);
                        temp = (data(j, i - 1) + data(j, i + 1)
                                + data(j - 1, i) + data(j + 1, i)) >> 2;
                        keep |= (copy - temp).toAbs() > limit;
                    }
                }
            }
            ++iterate;
        } while(keep && iterate < iterateMax);
    }

    template<class T>
    bool FDM<T>::onBoundary(size_t column, size_t row) {
        bool result = false;
        for(auto& boundary : boundaries) {
            if(boundary.type == Row)
                result |= row == boundary.index && boundary.from <= column && column <= boundary.to;
            else
                result |= column == boundary.index && boundary.from <= row && row <= boundary.to;
        }
        return result;
    }
}

#endif