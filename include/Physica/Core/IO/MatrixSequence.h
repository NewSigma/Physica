/*
 * Copyright 2022 WeiBo He.
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
#include <utility>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class ScalarType>
    class MatrixSequence {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector>;
        MatrixType mat;
        std::ifstream fin;
        uint64_t stepNum;
    public:
        MatrixSequence(size_t row, size_t column, std::ifstream fin_);
        MatrixSequence(const MatrixSequence&) = default;
        MatrixSequence(MatrixSequence&&) = default;
        ~MatrixSequence() = default;
        /* Operators */
        MatrixSequence& operator=(MatrixSequence obj) noexcept;
        /* Operations */
        bool step();
        /* Getters */
        [[nodiscard]] const MatrixType& getCurrent() const noexcept { return mat; }
        [[nodiscard]] uint64_t getStep() const noexcept { return stepNum; }
        /* Helpers */
        void swap(MatrixSequence& obj) noexcept;
    };

    template<class ScalarType>
    MatrixSequence<ScalarType>::MatrixSequence(size_t row, size_t column, std::ifstream fin_)
            : mat(row, column)
            , fin(std::move(fin_))
            , stepNum(0) {}

    template<class ScalarType>
    MatrixSequence<ScalarType>& MatrixSequence<ScalarType>::operator=(MatrixSequence obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * \returns true if read successfully
     */
    template<class ScalarType>
    bool MatrixSequence<ScalarType>::step() {
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        fin >> mat;
        ++stepNum;
        return bool(fin);
    }

    template<class ScalarType>
    void MatrixSequence<ScalarType>::swap(MatrixSequence& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mat.swap(obj.mat);
        fin.swap(obj.fin);
        std::swap(stepNum, obj.stepNum);
    }
}
