/*
 * Copyright 2021-2025 Weibo He.
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

#include "../LValueMatrix.h"

namespace Physica {
    namespace Internal {
        template<Matrix SourceType, Matrix TargetType, size_t Size>
        struct InverseImpl {
            static void run(const SourceType& source, TargetType& target) {
                const ssize_t order = source.getRow();
                const ssize_t order1 = ssize_t(order) - 1;
                SourceType copy = source;
                if constexpr (MatrixOption::isSameMajor<SourceType, TargetType>()) {
                    target.toUnitMatrix();
                    for (ssize_t i = 0; i < order1; ++i) {
                        ssize_t k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            ++k;
                            [[maybe_unused]] const bool isNotSingular = k < order;
                            assert(isNotSingular);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            target.majorSwap(k, i);
                        }

                        for (ssize_t j = i + 1; j < order; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            target.majorReduce(j, i, factor);
                        }
                    }

                    for (ssize_t i = order1; i > 0; --i) {
                        ssize_t k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            --k;
                            [[maybe_unused]] const bool isNotSingular = k < order;
                            assert(isNotSingular);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            target.majorSwap(k, i);
                        }

                        for (ssize_t j = 0; j < i; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            target.majorReduce(j, i, factor);
                        }
                    }
                    for (ssize_t i = 0; i < order; ++i)
                        target.majorMulScalar(i, reciprocal(copy(i, i)));
                }
                else {
                    auto temp = SourceType::unitMatrix(order);
                    for (ssize_t i = 0; i < order1; ++i) {
                        auto k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            ++k;
                            assert(k < order);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            temp.majorSwap(k, i);
                        }

                        for (ssize_t j = i + 1; j < order; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            temp.majorReduce(j, i, factor);
                        }
                    }

                    for (ssize_t i = order1; i > 0; --i) {
                        auto k = i;
                        while (copy.refFromMajorMinor(k, i).isZero()) {
                            --k;
                            assert(k < order);
                        }
                        if (k != i) {
                            copy.majorSwap(k, i);
                            target.majorSwap(k, i);
                        }

                        for (ssize_t j = 0; j < i; ++j) {
                            auto factor = copy.refFromMajorMinor(j, i) / copy.refFromMajorMinor(i, i);
                            copy.majorReduce(j, i, factor);
                            temp.majorReduce(j, i, factor);
                        }
                    }
                    for (ssize_t i = 0; i < order; ++i)
                        temp.majorMulScalar(i, reciprocal(copy(i, i)));
                    target = temp;
                }
            }
        };

        template<Matrix SourceType, Matrix TargetType>
        struct InverseImpl<SourceType, TargetType, 3> {
            static void run(const SourceType& source, TargetType& target) {
                const auto repDet = reciprocal(source.determinate());
                if constexpr (MatrixOption::isRowMatrix<TargetType>()) {
                    target(0, 0) = (source(1, 1) * source(2, 2) - source(1, 2) * source(2, 1)) * repDet;
                    target(0, 1) = (source(2, 1) * source(0, 2) - source(0, 1) * source(2, 2)) * repDet;
                    target(0, 2) = (source(0, 1) * source(1, 2) - source(1, 1) * source(0, 2)) * repDet;
                    target(1, 0) = (source(2, 0) * source(1, 2) - source(1, 0) * source(2, 2)) * repDet;
                    target(1, 1) = (source(0, 0) * source(2, 2) - source(2, 0) * source(0, 2)) * repDet;
                    target(1, 2) = (source(1, 0) * source(0, 2) - source(0, 0) * source(1, 2)) * repDet;
                    target(2, 0) = (source(1, 0) * source(2, 1) - source(2, 0) * source(1, 1)) * repDet;
                    target(2, 1) = (source(2, 0) * source(0, 1) - source(0, 0) * source(2, 1)) * repDet;
                    target(2, 2) = (source(0, 0) * source(1, 1) - source(1, 0) * source(0, 1)) * repDet;
                }
                else {
                    target(0, 0) = (source(1, 1) * source(2, 2) - source(1, 2) * source(2, 1)) * repDet;
                    target(1, 0) = (source(2, 0) * source(1, 2) - source(1, 0) * source(2, 2)) * repDet;
                    target(2, 0) = (source(1, 0) * source(2, 1) - source(2, 0) * source(1, 1)) * repDet;
                    target(0, 1) = (source(2, 1) * source(0, 2) - source(0, 1) * source(2, 2)) * repDet;
                    target(1, 1) = (source(0, 0) * source(2, 2) - source(2, 0) * source(0, 2)) * repDet;
                    target(2, 1) = (source(2, 0) * source(0, 1) - source(0, 0) * source(2, 1)) * repDet;
                    target(0, 2) = (source(0, 1) * source(1, 2) - source(1, 1) * source(0, 2)) * repDet;
                    target(1, 2) = (source(1, 0) * source(0, 2) - source(0, 0) * source(1, 2)) * repDet;
                    target(2, 2) = (source(0, 0) * source(1, 1) - source(1, 0) * source(0, 1)) * repDet;
                }
            }
        };
    }

    template<Matrix T>
    class InverseMatrix<T> : public RValueMatrix<InverseMatrix<T>> {
        const T& matrix;
    public:
        InverseMatrix(const LValueMatrix<T>& matrix_) : matrix(matrix_.getDerived()) {
            assert(matrix.getRow() == matrix.getCol());
        }
        template<Matrix M>
        void assign(M& target) const;
        /* Getters */
        [[nodiscard]] const T& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return matrix.getRow(); }
    };

    template<Matrix T>
    template<Matrix M>
    void InverseMatrix<T>::assign(M& target) const {
        constexpr size_t Size = T::RowAtCompile == Dynamic ? M::RowAtCompile : T::RowAtCompile;
        Internal::InverseImpl<T, M, Size>::run(matrix, target);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<InverseMatrix<T>> : public Traits<T> {};
}
