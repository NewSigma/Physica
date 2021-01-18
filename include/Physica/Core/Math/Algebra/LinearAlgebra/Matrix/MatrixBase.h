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
namespace Physica::Core {
    namespace Intenal {
        template<class T>
        class MatrixTrait;
        /**
         * The \class MatrixBase provide algorithms that a matrix should support, but the data is provided
         * by \tparam MatrixImpl.
         * 
         * \tparam MatrixImpl
         * A class that contains data structure for a matrix.
         */
        template<class MatrixImpl>
        class MatrixBase {
        public:
            using ScalarType = typename MatrixTrait<MatrixImpl>::ScalarType;
        public:
            ScalarType determinate() const;
        private:
            MatrixImpl& getImpl() { return static_cast<MatrixImpl&>(*this); }
            const MatrixImpl& getImpl() const { return static_cast<MatrixImpl&>(*this); }
        };
    }
}