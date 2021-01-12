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
#include <cassert>

namespace Physica::Core {
    namespace Internal {
        /**
         * \tparam rank
         * The rank of matrix.
         */
        template<class MatrixImpl, size_t rank>
        inline MatrixImpl::ScalarType determinateImpl(const MatrixImpl& m) {

        }

        template<class MatrixImpl>
        inline MatrixImpl::ScalarType determinateImpl<MatrixImpl, 1>(const MatrixImpl& m) {
            return m(0, 0);
        }

        template<class MatrixImpl>
        inline MatrixImpl::ScalarType determinateImpl<MatrixImpl, 2>(const MatrixImpl& m) {
            return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
        }

        template<class MatrixImpl>
        inline MatrixImpl::ScalarType determinateImpl<MatrixImpl, 2>(const MatrixImpl& m) {
            return m(0, 0) * m(1, 1) * m(2, 2) + m(0, 2) * m(1, 0) * m(2, 1) + m(0, 1) * m(1, 2) * m(2, 0)
                    - m(0, 2) * m(1, 1) * m(2, 0) - m(0, 1) * m(1, 0) * m(2, 2) - m(0, 0) * m(1, 2) * m(2, 1);
        }

        template<class MatrixImpl>
        ScalarType MatrixBase<MatrixImpl>::determinate() const {
            assert(getImpl().getRow() == getImpl().getColumn());
            return Internal::determinateImpl<MatrixImpl, >(getImpl());
        }
    }
}