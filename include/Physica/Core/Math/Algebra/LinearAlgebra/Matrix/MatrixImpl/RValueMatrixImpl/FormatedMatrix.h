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
    template<Matrix M>
    class FormatedMatrix {
        using This = FormatedMatrix<M>;
        using ScalarType = M::ScalarType;

        const M& data;
        std::string matPrefix;
        std::string matSuffix;
        std::string rowPrefix;
        std::string rowSuffix;
        std::string rowSeparator;
        std::string separator;
    public:
        FormatedMatrix(const M& data_);
        FormatedMatrix(const This&) = delete;
        FormatedMatrix(This&&) noexcept = delete;
        ~FormatedMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        This& toFormatMMA();
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
        This& setMatPrefix(std::string matPrefix_);
        This& setMatSuffix(std::string matSuffix_);
        This& setRowPrefix(std::string rowPrefix_);
        This& setRowSuffix(std::string rowSuffix_);
        This& setRowSeparator(std::string rowSeparator_);
        This& setSeparator(std::string separator_);
    private:
        [[nodiscard]] bool isVectorLike() const noexcept { return (getRow() == 1) || (getCol() == 1); }
    };

    template<Matrix M>
    FormatedMatrix<M>::FormatedMatrix(const M& data_)
            : data(data_)
            , matSuffix("\n")
            , rowSeparator("\n")
            , separator(" ") {
        if (isVectorLike()) {
            setMatPrefix("(");
            setMatSuffix(")");
            setSeparator(", ");
        }
    }

    template<Matrix M>
    auto FormatedMatrix<M>::toFormatMMA() -> This& {
        setMatPrefix("{");
        setMatSuffix("}");
        setRowPrefix("{");
        setRowSuffix("}");
        setRowSeparator(",");
        setSeparator(",");
        return *this;
    }

    template<Matrix M>
    auto FormatedMatrix<M>::setMatPrefix(std::string matPrefix_) -> This& {
        matPrefix = std::move(matPrefix_);
        return *this;
    }

    template<Matrix M>
    auto FormatedMatrix<M>::setMatSuffix(std::string matSuffix_) -> This& {
        matSuffix = std::move(matSuffix_);
        return *this;
    }

    template<Matrix M>
    auto FormatedMatrix<M>::setRowPrefix(std::string rowPrefix_) -> This& {
        rowPrefix = std::move(rowPrefix_);
        return *this;
    }

    template<Matrix M>
    auto FormatedMatrix<M>::setRowSuffix(std::string rowSuffix_) -> This& {
        rowSuffix = std::move(rowSuffix_);
        return *this;
    }

    template<Matrix M>
    auto FormatedMatrix<M>::setRowSeparator(std::string rowSeparator_) -> This& {
        rowSeparator = std::move(rowSeparator_);
        return *this;
    }

    template<Matrix M>
    auto FormatedMatrix<M>::setSeparator(std::string separator_) -> This& {
        separator = std::move(separator_);
        return *this;
    }

    template<Matrix U>
    std::ostream& operator<<(std::ostream& os, const FormatedMatrix<U>& m) {
        return os << std::format("{}", m);
    }
}

namespace std {
    template<Physica::Matrix M>
    struct formatter<Physica::FormatedMatrix<M>, char> {
        constexpr auto parse(std::format_parse_context& ctx);
        static auto format(const Physica::FormatedMatrix<M>& obj, std::format_context& ctx) {
            const size_t col = obj.getCol();
            const size_t row = obj.getRow();
            const auto width = [&obj, row, col]() noexcept {
                size_t result = 0;
                for (size_t c = 0; c < col; ++c)
                    for (size_t r = 0; r < row; ++r)
                        result = std::max(result, std::formatted_size("{}", obj.calc(r, c).real()));
                return static_cast<std::streamsize>(result);
            }();

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
                ss << obj.getRowSuffix();

                const bool isLastRow = r + 1 == row;
                if (!isLastRow)
                    ss << obj.getRowSeparator();
            }
            ss << obj.getMatSuffix();
            return std::format_to(ctx.out(), "{}", ss.str());
        }
    };
}
