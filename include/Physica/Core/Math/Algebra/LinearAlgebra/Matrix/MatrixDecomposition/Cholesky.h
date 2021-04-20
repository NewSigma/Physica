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

        template<bool isVectorMatrix, class MatrixIterator>
        struct GetElementIteratorImpl {
            inline static auto run(const MatrixIterator& matrixIte) {
                return (*matrixIte).begin();
            }
        };

        template<class MatrixIterator>
        struct GetElementIteratorImpl<false, MatrixIterator> {
            inline static auto run(const MatrixIterator& matrixIte) {
                return matrixIte;
            }
        };

        template<bool isVectorMatrix, class MatrixIterator>
        struct GetConstElementIteratorImpl {
            inline static auto run(const MatrixIterator& matrixIte) {
                return (*matrixIte).cbegin();
            }
        };

        template<class MatrixIterator>
        struct GetConstElementIteratorImpl<false, MatrixIterator> {
            inline static auto run(const MatrixIterator& matrixIte) {
                return matrixIte;
            }
        };

        template<bool isVectorMatrix, class MatrixIterator, class ElementIterator>
        struct UpdateIteratorImpl {
            inline static void run(MatrixIterator& matrixIte, ElementIterator& elementIte) {
                ++matrixIte;
                elementIte = (*matrixIte).begin();
            }
        };

        template<class MatrixIterator, class ElementIterator>
        struct UpdateIteratorImpl<false, MatrixIterator, ElementIterator> {
            inline static void run([[maybe_unused]] MatrixIterator& matrixIte, ElementIterator& elementIte) {
                ++elementIte;
            }
        };

        template<bool isVectorMatrix, class MatrixIterator, class ElementIterator>
        struct UpdateConstIteratorImpl {
            inline static void run(MatrixIterator& matrixIte, ElementIterator& elementIte) {
                ++matrixIte;
                elementIte = (*matrixIte).cbegin();
            }
        };

        template<class MatrixIterator, class ElementIterator>
        struct UpdateConstIteratorImpl<false, MatrixIterator, ElementIterator> {
            inline static void run([[maybe_unused]] MatrixIterator& matrixIte, ElementIterator& elementIte) {
                ++elementIte;
            }
        };
    }

    template<class Matrix>
    class Cholesky {
        const Matrix& matrix;
    public:
        explicit Cholesky(const Matrix& matrix_) : matrix(matrix_) {}
        ~Cholesky() = default;
        /* Getters */
        [[nodiscard]] const Matrix& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] size_t getOrder() const noexcept { return matrix.getOrder(); }
        /* Static members */
        template<class OtherMatrix, class MatrixIterator>
        inline static auto getElementIterator(const MatrixIterator& matrixIte) {
            constexpr int isVectorMatrix = DenseMatrixType::isVectorMatrix(Internal::Traits<OtherMatrix>::MatrixType);
            return Internal::GetElementIteratorImpl<isVectorMatrix, MatrixIterator>::run(matrixIte);
        }

        template<class OtherMatrix, class MatrixIterator>
        inline static auto getConstElementIterator(const MatrixIterator& matrixIte) {
            constexpr int isVectorMatrix = DenseMatrixType::isVectorMatrix(Internal::Traits<OtherMatrix>::MatrixType);
            return Internal::GetConstElementIteratorImpl<isVectorMatrix, MatrixIterator>::run(matrixIte);
        }

        template<class OtherMatrix, class MatrixIterator, class ElementIterator>
        inline static void updateConstIterator(MatrixIterator& matrixIte, ElementIterator& elementIte) {
            constexpr int isVectorMatrix = DenseMatrixType::isVectorMatrix(Internal::Traits<OtherMatrix>::MatrixType);
            Internal::UpdateConstIteratorImpl<isVectorMatrix, MatrixIterator, ElementIterator>::run(matrixIte, elementIte);
        }

        template<class OtherMatrix, class MatrixIterator, class ElementIterator>
        inline static void updateIterator(MatrixIterator& matrixIte, ElementIterator& elementIte) {
            constexpr int isVectorMatrix = DenseMatrixType::isVectorMatrix(Internal::Traits<OtherMatrix>::MatrixType);
            Internal::UpdateIteratorImpl<isVectorMatrix, MatrixIterator, ElementIterator>::run(matrixIte, elementIte);
        }
    };
}
