/*
 * Copyright 2024 Weibo He.
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
     * \class FormatedMatrix convert a matrix to text, either readable to human, or other softwares.
     */
    template<class MatrixType>
    class FormatedMatrix {
        using ScalarType = typename MatrixType::ScalarType;

        const RValueMatrix<MatrixType>& data;
        std::string matPrefix;
        std::string matSuffix;
        std::string rowPrefix;
        std::string rowSuffix;
        std::string rowSeparator;
        std::string separator;
    public:
        FormatedMatrix(const RValueMatrix<MatrixType>& data_);
        /* Operators */
        template<class T>
        friend std::ostream& operator<<(std::ostream& os, const FormatedMatrix<T>& m);
        /* Operations */
        FormatedMatrix& toFormatMMA();
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return data.calc(row, col); }
        [[nodiscard]] size_t getRow() const noexcept { return data.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return data.getCol(); }
        /* Setters */
        FormatedMatrix& setMatPrefix(std::string matPrefix_);
        FormatedMatrix& setMatSuffix(std::string matSuffix_);
        FormatedMatrix& setRowPrefix(std::string rowPrefix_);
        FormatedMatrix& setRowSuffix(std::string rowSuffix_);
        FormatedMatrix& setRowSeparator(std::string rowSeparator_);
        FormatedMatrix& setSeparator(std::string separator_);
    };

    template<class MatrixType>
    FormatedMatrix<MatrixType>::FormatedMatrix(const RValueMatrix<MatrixType>& data_)
            : data(data_)
            , matPrefix("")
            , matSuffix("")
            , rowPrefix("")
            , rowSuffix("")
            , rowSeparator("")
            , separator(" ") {}

    template<class MatrixType>
    std::ostream& operator<<(std::ostream& os, const FormatedMatrix<MatrixType>& m) {
        const size_t col = m.getCol();
        const size_t row = m.getRow();
        size_t width = 0;
        /* Get max width */ {
            for (size_t c = 0; c < col; ++c) {
                for (size_t r = 0; r < row; ++ r) {
                    std::stringstream stream{};
                    stream.copyfmt(os);
                    stream << m.calc(r, c).real();
                    width = std::max(width, stream.str().length());
                }
            }
        }
        /* Output */ {
            os << m.matPrefix;
            for (size_t r = 0; r < row; ++r) {
                os << m.rowPrefix;
                for (size_t c = 0; c < col; ++c) {
                    os.width(width);
                    os << m.calc(r, c);
                    const bool isLastElem = c + 1 == col;
                    if (!isLastElem)
                        os << m.separator;
                }
                os << m.rowSuffix << '\n';

                const bool isLastRow = r + 1 == row;
                if (!isLastRow)
                    os << m.rowSeparator;
            }
            os << m.matSuffix;
        }
        return os;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::toFormatMMA() {
        setMatPrefix("{");
        setMatSuffix("}");
        setRowPrefix("{");
        setRowSuffix("}");
        setRowSeparator(",");
        setSeparator(",");
        return *this;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::setMatPrefix(std::string matPrefix_) {
        matPrefix = std::move(matPrefix_);
        return *this;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::setMatSuffix(std::string matSuffix_) {
        matSuffix = std::move(matSuffix_);
        return *this;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::setRowPrefix(std::string rowPrefix_) {
        rowPrefix = std::move(rowPrefix_);
        return *this;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::setRowSuffix(std::string rowSuffix_) {
        rowSuffix = std::move(rowSuffix_);
        return *this;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::setRowSeparator(std::string rowSeparator_) {
        rowSeparator = std::move(rowSeparator_);
        return *this;
    }

    template<class MatrixType>
    FormatedMatrix<MatrixType>& FormatedMatrix<MatrixType>::setSeparator(std::string separator_) {
        separator = std::move(separator_);
        return *this;
    }
}
