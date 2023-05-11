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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"

namespace Physica::Core {
    template<class MatrixType, bool isLValueMatrix> class DiagVector;

    namespace Internal {
        template<class MatrixType, bool isLValueMatrix>
        class Traits<DiagVector<MatrixType, isLValueMatrix>> {
            static_assert(std::is_convertible<MatrixType&, Core::LValueMatrix<MatrixType>&>::value == isLValueMatrix, "[Error]: Invalid LValueMatrix");
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile > MatrixType::ColumnAtCompile ? MatrixType::RowAtCompile : MatrixType::ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = MatrixType::MaxRowAtCompile > MatrixType::MaxColumnAtCompile ? MatrixType::MaxRowAtCompile : MatrixType::MaxColumnAtCompile;
        };
    }

    template<class MatrixType, bool isLValueMatrix>
    class DiagVector : public RValueVector<DiagVector<MatrixType, isLValueMatrix>> {
        const MatrixType& mat;
        using Base = RValueVector<DiagVector<MatrixType, isLValueMatrix>>;
    public:
        using typename Base::ScalarType;
    public:
        explicit DiagVector(const MatrixType& mat_) : mat(mat_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return mat.calc(index, index); }
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
    };

    template<class MatrixType>
    class DiagVector<MatrixType, true> : public LValueVector<DiagVector<MatrixType, true>> {
        using Base = LValueVector<DiagVector<MatrixType, true>>;
    public:
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
    public:
        explicit DiagVector(MatrixType& mat_) : mat(mat_) {}
        /* Operations */
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] ScalarType& operator[](size_t index) { return mat(index, index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return mat(index, index); }
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
    };
}
