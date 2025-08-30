/*
 * Copyright 2024-2025 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    /**
     * \class FormatedMatrix convert a matrix to text, either readable to human, or other softwares.
     */
    template<Matrix T>
    class FormatedMatrix {
        using ScalarType = T::ScalarType;

        const T& data;
        std::string matPrefix;
        std::string matSuffix;
        std::string rowPrefix;
        std::string rowSuffix;
        std::string rowSeparator;
        std::string separator;
    public:
        FormatedMatrix(const T& data_);
        /* Operations */
        FormatedMatrix& toFormatMMA();
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return data.calc(row, col); }
        [[nodiscard]] size_t getRow() const noexcept { return data.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return data.getCol(); }
        [[nodiscard]] const std::string& getMatPrefix() const noexcept { return matPrefix; }
        [[nodiscard]] const std::string& getMatSuffix() const noexcept { return matSuffix; }
        [[nodiscard]] const std::string& getRowPrefix() const noexcept { return rowPrefix; }
        [[nodiscard]] const std::string& getRowSuffix() const noexcept { return rowSuffix; }
        [[nodiscard]] const std::string& getRowSeparator() const noexcept { return rowSeparator; }
        [[nodiscard]] const std::string& getSeparator() const noexcept { return separator; }
        /* Setters */
        FormatedMatrix& setMatPrefix(std::string matPrefix_);
        FormatedMatrix& setMatSuffix(std::string matSuffix_);
        FormatedMatrix& setRowPrefix(std::string rowPrefix_);
        FormatedMatrix& setRowSuffix(std::string rowSuffix_);
        FormatedMatrix& setRowSeparator(std::string rowSeparator_);
        FormatedMatrix& setSeparator(std::string separator_);
    };

    template<Matrix T>
    FormatedMatrix<T>::FormatedMatrix(const T& data_)
            : data(data_)
            , matPrefix("")
            , matSuffix("")
            , rowPrefix("")
            , rowSuffix("")
            , rowSeparator("")
            , separator(" ") {
        bool isVectorLike = (getRow() == 1) || (getCol() == 1);
        if (isVectorLike) {
            setMatPrefix("(");
            setMatSuffix(")");
            setSeparator(", ");
        }
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::toFormatMMA() {
        setMatPrefix("{");
        setMatSuffix("}");
        setRowPrefix("{");
        setRowSuffix("}");
        setRowSeparator(",");
        setSeparator(",");
        return *this;
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::setMatPrefix(std::string matPrefix_) {
        matPrefix = std::move(matPrefix_);
        return *this;
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::setMatSuffix(std::string matSuffix_) {
        matSuffix = std::move(matSuffix_);
        return *this;
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::setRowPrefix(std::string rowPrefix_) {
        rowPrefix = std::move(rowPrefix_);
        return *this;
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::setRowSuffix(std::string rowSuffix_) {
        rowSuffix = std::move(rowSuffix_);
        return *this;
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::setRowSeparator(std::string rowSeparator_) {
        rowSeparator = std::move(rowSeparator_);
        return *this;
    }

    template<Matrix T>
    FormatedMatrix<T>& FormatedMatrix<T>::setSeparator(std::string separator_) {
        separator = std::move(separator_);
        return *this;
    }

    template<Matrix U>
    std::ostream& operator<<(std::ostream& os, const FormatedMatrix<U>& m) {
        return os << std::format("{}", m);
    }
}

namespace std {
    template<Physica::Matrix T>
    struct formatter<Physica::FormatedMatrix<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::FormatedMatrix<T>& obj, std::format_context& ctx) const {
            const size_t col = obj.getCol();
            const size_t row = obj.getRow();

            size_t width = 0;
            for (size_t c = 0; c < col; ++c)
                for (size_t r = 0; r < row; ++ r)
                    width = std::max(width, std::formatted_size("{}", obj.calc(r, c).real()));

            std::stringstream ss{};
            ss << obj.getMatPrefix();
            for (size_t r = 0; r < row; ++r) {
                ss << obj.getRowPrefix();
                for (size_t c = 0; c < col; ++c) {
                    ss.width(width);
                    ss << obj.calc(r, c);
                    const bool isLastElem = c + 1 == col;
                    if (!isLastElem)
                        ss << obj.getSeparator();
                }
                ss << obj.getRowSuffix() << '\n';

                const bool isLastRow = r + 1 == row;
                if (!isLastRow)
                    ss << obj.getRowSeparator();
            }
            ss << obj.getMatSuffix();
            return std::format_to(ctx.out(), "{}", ss.str());
        }
    };
}
