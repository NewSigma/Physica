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

namespace Physica::Core {
    /**
     * Decomposite a symmetrical, positive matrix A into LL^T.
     * If target matrix is a column matrix, return lower triangular matrix L, if row matrix, return upper triangular matrix L^T
     */
    template<Matrix T>
    class Cholesky : public RValueMatrix<Cholesky<T>> {
        using Base = RValueMatrix<Cholesky<T>>;
        const T& source;
    public:
        explicit Cholesky(const T& source_);
        ~Cholesky() = default;
        /* Operations */
        template<Matrix M>
        void assignTo(LValueMatrix<M>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return source.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return source.getRow(); }
    };

    template<Matrix T>
    Cholesky<T>::Cholesky(const T& source_) : source(source_) {
        assert(source.getRow() == source.getCol());
    }
    /**
     * Implemented the square method
     */
    template<Matrix T>
    template<Matrix M>
    void Cholesky<T>::assignTo(LValueMatrix<M>& target) const {
        using ResultType = M;
        using ScalarType = typename ResultType::ScalarType;
        const size_t order = source.getRow();
        /* Handle first vector */ {
            const auto diag = sqrt(source(0, 0));
            target(0, 0) = diag;
            const ScalarType inv_diag = reciprocal(diag);
            for (size_t minor = 1; minor < order; ++minor)
                target.refFromMajorMinor(0, minor) = source.refFromMajorMinor(0, minor) * inv_diag;
        }
        /* Handle other vectors */ {
            for (size_t major = 1; major < order; ++major) {
                size_t minor = 0;
                for (; minor < major; ++minor)
                    target.refFromMajorMinor(major, minor) = ScalarType(0);

                ScalarType diag = source(major, major);
                /* major == minor */ {
                    for (size_t k = 0; k < major; ++k)
                        diag -= square(target.refFromMajorMinor(k, major));
                    diag = sqrt(diag);
                    target(major, major) = diag;
                    ++minor;
                }
                const ScalarType inv_diag = reciprocal(diag);

                for (; minor < order; ++minor) {
                    ScalarType temp = source.refFromMajorMinor(major, minor);
                    for (size_t k = 0; k < major; ++k)
                        temp -= target.refFromMajorMinor(k, major) * target.refFromMajorMinor(k, minor);
                    target.refFromMajorMinor(major, minor) = temp * inv_diag;
                }
            }
        }
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<Cholesky<T>> {
    public:
        using ScalarType = typename T::ScalarType;
        constexpr static int Option = T::MatrixOption;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
