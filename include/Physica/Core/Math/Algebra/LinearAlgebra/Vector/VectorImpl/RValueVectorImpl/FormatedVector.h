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

#include "../RValueVector.h"

namespace Physica {
    template<class VectorType> class FormatedVector;

    template<Vector T>
    class FormatedVector<T> {
        using This = FormatedVector<T>;

        const T& data;
        std::string prefix;
        std::string suffix;
        std::string separator;
    public:
        FormatedVector(const T& data_);
        FormatedVector(const This&) = delete;
        FormatedVector(This&&) noexcept = delete;
        ~FormatedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        This& toFormatMMA();
        /* Getters */
        [[nodiscard]] const T& getData() const noexcept { return data; }
        [[nodiscard]] const std::string& getPrefix() const noexcept { return prefix; }
        [[nodiscard]] const std::string& getSuffix() const noexcept { return suffix; }
        [[nodiscard]] const std::string& getSeparator() const noexcept { return separator; }
        /* Setters */
        This& setPrefix(std::string prefix_);
        This& setSuffix(std::string suffix_);
        This& setSeparator(std::string separator_);
    };

    template<Vector T>
    FormatedVector<T>::FormatedVector(const T& data_)
            : data(data_)
            , prefix("{")
            , suffix("}")
            , separator(", ") {}

    template<Vector T>
    auto FormatedVector<T>::toFormatMMA() -> This& {
        setPrefix("{");
        setSuffix("}");
        setSeparator(",");
        return *this;
    }

    template<Vector T>
    auto FormatedVector<T>::setPrefix(std::string prefix_) -> This& {
        prefix = std::move(prefix_);
        return *this;
    }

    template<Vector T>
    auto FormatedVector<T>::setSuffix(std::string suffix_) -> This& {
        suffix = std::move(suffix_);
        return *this;
    }

    template<Vector T>
    auto FormatedVector<T>::setSeparator(std::string separator_) -> This& {
        separator = std::move(separator_);
        return *this;
    }

    template<Vector U>
    std::ostream& operator<<(std::ostream& os, const FormatedVector<U>& v) {
        return os << std::format("{}", v);
    }
}

namespace std {
    template<Physica::Vector T>
    struct formatter<Physica::FormatedVector<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); };
        static auto format(const Physica::FormatedVector<T>& obj, auto& ctx) {
            const auto& data = obj.getData();
            const auto& sep = obj.getSeparator();
            size_t length = data.getLength();
            std::stringstream ss{};
            ss << obj.getPrefix();
            if (length > 0) {
                --length;
                for (size_t i = 0; i < length; ++i)
                    ss << data.calc(i) << sep;
                ss << data.calc(length);
            }
            ss << obj.getSuffix();
            return std::format_to(ctx.out(), "{}", ss.str());
        }
    };
}
