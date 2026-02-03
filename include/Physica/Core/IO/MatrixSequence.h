/*
 * Copyright 2022-2025 Weibo He.
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

namespace Physica {
    template<Scalar T>
    class MatrixSequence {
        using This = MatrixSequence<T>;
        using MatrixType = DenseMatrix<T, MatrixMajor::Row>;
        MatrixType mat;
        std::ifstream fin;
        uint64_t stepNum;
    public:
        MatrixSequence(size_t row, size_t col, std::ifstream fin_);
        MatrixSequence(const This&) = default;
        MatrixSequence(This&&) = default;
        ~MatrixSequence() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        bool step();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const MatrixType& getCurrent() const noexcept { return mat; }
        [[nodiscard]] uint64_t getStep() const noexcept { return stepNum; }
    };

    template<Scalar T>
    MatrixSequence<T>::MatrixSequence(size_t row, size_t col, std::ifstream fin_)
            : mat(row, col)
            , fin(std::move(fin_))
            , stepNum(0) {}
    /**
     * \returns true if read successfully
     */
    template<Scalar T>
    bool MatrixSequence<T>::step() {
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        fin >> mat;
        ++stepNum;
        return bool(fin);
    }

    template<Scalar T>
    void MatrixSequence<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mat.swap(obj.mat);
        fin.swap(obj.fin);
        std::swap(stepNum, obj.stepNum);
    }
}
