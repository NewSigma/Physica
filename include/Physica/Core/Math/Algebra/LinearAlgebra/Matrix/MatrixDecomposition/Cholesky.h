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
    namespace Internal {
        template<class T>
        class Traits;
    }

    template<class Matrix>
    class Cholesky {
        const Matrix& matrix;
    public:
        explicit Cholesky(const Matrix& matrix);
        ~Cholesky() = default;
        /* Operations */
        //A line may be a row or a column, depending on row major or column major
        void calculateLine(size_t lineNumber);
        /* Getters */
        [[nodiscard]] const Matrix& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] size_t getOrder() const noexcept { return matrix.getOrder(); }
    };
}
