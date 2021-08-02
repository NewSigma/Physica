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
    template<class MatrixType> class Cholesky;

    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<Cholesky<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static int MatrixOption = MatrixType::MatrixOption;
            constexpr static size_t RowAtCompile = MatrixType::RowAtCompile;
            constexpr static size_t ColumnAtCompile = MatrixType::ColumnAtCompile;
            constexpr static size_t MaxRowAtCompile = MatrixType::MaxRowAtCompile;
            constexpr static size_t MaxColumnAtCompile = MatrixType::MaxColumnAtCompile;
            constexpr static size_t SizeAtCompile = MatrixType::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = MatrixType::MaxSizeAtCompile;
        };
    }

    template<class MatrixType>
    class Cholesky : public RValueMatrix<Cholesky<MatrixType>> {
        using Base = RValueMatrix<Cholesky<MatrixType>>;
        const MatrixType& matrix;
    public:
        explicit Cholesky(const MatrixType& matrix_);
        ~Cholesky() = default;
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& mat) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return matrix.getRow(); }
    };

    template<class MatrixType>
    Cholesky<MatrixType>::Cholesky(const MatrixType& matrix_) : matrix(matrix_) {
        assert(matrix.getRow() == matrix.getColumn());
    }
    /**
     * Implemented the square method
     */
    template<class MatrixType>
    template<class OtherMatrix>
    void Cholesky<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& mat) const {
        using ResultType = OtherMatrix;
        using ScalarType = typename ResultType::ScalarType;
        const size_t order = matrix.getRow();
        auto matrixIte = mat.getDerived().begin();
        auto constMatrixIte = matrix.cbegin();

        auto elementIte = mat.getDerived().ebegin(matrixIte);
        auto constElementIte = matrix.cebegin(constMatrixIte);
        /* Handle first vector */ {
            const auto diag = sqrt(*constElementIte);
            *elementIte = diag;
            for (size_t i = 1; i < order; ++i) {
                ++elementIte;
                ++constElementIte;
                *elementIte = *constElementIte / diag;
            }
        }
        /* Handle other vectors */ {
            for (size_t i = 1; i < order; ++i) {
                ResultType::updateIterator(matrixIte, elementIte);
                MatrixType::updateIterator(constMatrixIte, constElementIte);
                size_t j;
                for (j = 0; j < i; ++j) {
                    *elementIte = 0;
                    ++elementIte;
                    ++constElementIte;
                }

                ScalarType diag(*constElementIte);
                /* j == i */ {
                    for (size_t k = 0; k < i; ++k) {
                        if constexpr (DenseMatrixOption::isColumnMatrix<ResultType>())
                            diag -= square(mat(i, k));
                        else
                            diag -= square(mat(k, i));
                    }
                    diag = sqrt(diag);
                    *elementIte = diag;
                    ++j;
                }

                for (; j < order; ++j) {
                    ++elementIte;
                    ++constElementIte;
                    ScalarType temp(*constElementIte);
                    for (size_t k = 0; k < j; ++k) {
                        if constexpr (DenseMatrixOption::isColumnMatrix<ResultType>())
                            diag -= mat(i, k) * mat(j, k);
                        else
                            diag -= mat(k, i) * mat(k, j);
                    }
                    *elementIte = temp / diag;
                }
            }
        }
    }
}
