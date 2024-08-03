/*
 * Copyright 2021-2024 Weibo He.
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

#include <cassert>
#include "MatrixImpl/LValueMatrix.h"

namespace Physica::Core {
    namespace Internal {
        template<class SourceType, class TargetType, size_t Size>
        struct InverseImpl {
            static void run(const LValueMatrix<SourceType>& source, LValueMatrix<TargetType>& target) {
                const size_t order = source.getRow();
                const size_t order_1 = order - 1;
                SourceType copy = source.getDerived();
                if constexpr (MatrixOption::isSameMajor<SourceType, TargetType>()) {
                    target.getDerived().toUnitMatrix();
                    for (size_t i = 0; i < order_1; ++i) {
                        size_t k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            ++k;
                            [[maybe_unused]] const bool isNotSingular = k < order;
                            assert(isNotSingular);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            target.getDerived().majorSwap(k, i);
                        }

                        for (size_t j = i + 1; j < order; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            target.getDerived().majorReduce(j, i, factor);
                        }
                    }

                    for (size_t i = order_1; i > 0; --i) {
                        size_t k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            --k;
                            [[maybe_unused]] const bool isNotSingular = k < order;
                            assert(isNotSingular);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            target.getDerived().majorSwap(k, i);
                        }

                        for (size_t j = 0; j < i; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            target.getDerived().majorReduce(j, i, factor);
                        }
                    }
                    for (size_t i = 0; i < order; ++i)
                        target.getDerived().majorMulScalar(i, reciprocal(copy(i, i)));
                }
                else {
                    auto temp = SourceType::unitMatrix(order);
                    for (size_t i = 0; i < order_1; ++i) {
                        size_t k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            ++k;
                            assert(k < order);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            temp.majorSwap(k, i);
                        }

                        for (size_t j = i + 1; j < order; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            temp.majorReduce(j, i, factor);
                        }
                    }

                    for (size_t i = order_1; i > 0; --i) {
                        size_t k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            --k;
                            assert(k < order);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            target.getDerived().majorSwap(k, i);
                        }

                        for (size_t j = 0; j < i; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            temp.majorReduce(j, i, factor);
                        }
                    }
                    for (size_t i = 0; i < order; ++i)
                        temp.majorMulScalar(i, reciprocal(copy(i, i)));
                    target.getDerived() = temp;
                }
            }
        };

        template<class SourceType, class TargetType>
        struct InverseImpl<SourceType, TargetType, 3> {
            static void run(const LValueMatrix<SourceType>& source, LValueMatrix<TargetType>& target) {
                using ScalarType = typename SourceType::ScalarType;
                const ScalarType repDet = reciprocal(source.getDerived().determinate());
                if constexpr (TargetType::isRowMatrix) {
                    target(0, 0) = (source(1, 1) * source(2, 2) - source(1, 2) * source(2, 1)) * repDet;
                    target(0, 1) = -(source(0, 1) * source(2, 2) - source(2, 1) * source(0, 2)) * repDet;
                    target(0, 2) = (source(0, 1) * source(1, 2) - source(1, 1) * source(0, 2)) * repDet;
                    target(1, 0) = -(source(1, 0) * source(2, 2) - source(2, 0) * source(1, 2)) * repDet;
                    target(1, 1) = (source(0, 0) * source(2, 2) - source(2, 0) * source(0, 2)) * repDet;
                    target(1, 2) = -(source(0, 0) * source(1, 2) - source(1, 0) * source(0, 2)) * repDet;
                    target(2, 0) = (source(1, 0) * source(2, 1) - source(2, 0) * source(1, 1)) * repDet;
                    target(2, 1) = -(source(0, 0) * source(2, 1) - source(2, 0) * source(0, 1)) * repDet;
                    target(2, 2) = (source(0, 0) * source(1, 1) - source(1, 0) * source(0, 1)) * repDet;
                }
                else {
                    target(0, 0) = (source(1, 1) * source(2, 2) - source(1, 2) * source(2, 1)) * repDet;
                    target(1, 0) = -(source(1, 0) * source(2, 2) - source(2, 0) * source(1, 2)) * repDet;
                    target(2, 0) = (source(1, 0) * source(2, 1) - source(2, 0) * source(1, 1)) * repDet;
                    target(0, 1) = -(source(0, 1) * source(2, 2) - source(2, 1) * source(0, 2)) * repDet;
                    target(1, 1) = (source(0, 0) * source(2, 2) - source(2, 0) * source(0, 2)) * repDet;
                    target(2, 1) = -(source(0, 0) * source(2, 1) - source(2, 0) * source(0, 1)) * repDet;
                    target(0, 2) = (source(0, 1) * source(1, 2) - source(1, 1) * source(0, 2)) * repDet;
                    target(1, 2) = -(source(0, 0) * source(1, 2) - source(1, 0) * source(0, 2)) * repDet;
                    target(2, 2) = (source(0, 0) * source(1, 1) - source(1, 0) * source(0, 1)) * repDet;
                }

                constexpr bool isContinuous = std::is_base_of<ContinuousMatrix<TargetType>, TargetType>::value;
                if constexpr (isContinuous && ScalarType::isReverseDiff)
                    target.getDerived().makeContinuous();
            }
        };
    }

    template<class MatrixType>
    class InverseMatrix : public RValueMatrix<InverseMatrix<MatrixType>> {
        const MatrixType& matrix;
    public:
        InverseMatrix(const LValueMatrix<MatrixType>& matrix_) : matrix(matrix_.getDerived()) {
            assert(matrix.getRow() == matrix.getColumn());
        }
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return matrix.getRow(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void InverseMatrix<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        constexpr size_t Size = MatrixType::RowAtCompile == Dynamic ? OtherMatrix::RowAtCompile : MatrixType::RowAtCompile;
        Internal::InverseImpl<MatrixType, OtherMatrix, Size>::run(matrix, target.getDerived());
    }
}

namespace Physica {
    template<class MatrixType>
    class Traits<Core::InverseMatrix<MatrixType>> : public Traits<MatrixType> {};
}
