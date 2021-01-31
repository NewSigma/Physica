/*
 * Copyright 2021 WeiBo He.
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

#include <iosfwd>
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    /**
     * AbstractMatrix is father class of DenseMatrix and SparseMatrix.
     */
    template<class Derived>
    class AbstractMatrix : public virtual Utils::CRTPBase<Derived> {
    public:
        std::ostream operator<<(std::ostream& os, const AbstractMatrix<Derived>& m);
    };

    std::ostream AbstractMatrix::operator<<(std::ostream& os, const AbstractMatrix<Derived>& m) {
        const auto row = m.getDerived().getRow();
        const auto column = m.getDerived().getColumn();
        //10 is the max precision of double.
        os << std::setprecision(10);
        for(size_t i = 0; i < row; ++i) {
            for(size_t j = 0; j < column; ++j)
                os << m.getDerived().operator()(row, column) << '\t';
            os << '\n';
        }
        //6 is the default precision.
        return os << std::setprecision(6);
    }
}
