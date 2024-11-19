/*
 * Copyright 2022 Weibo He.
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
    template<Vector T>
    class FormatedVector<T> {
        const T& data;
        std::string prefix;
        std::string suffix;
        std::string separator;
    public:
        FormatedVector(const T& data_);
        /* Operators */
        template<Vector U>
        friend std::ostream& operator<<(std::ostream& os, const FormatedVector<U>& v);
        /* Operations */
        FormatedVector& toFormatMMA();
        /* Setters */
        FormatedVector& setPrefix(std::string prefix_);
        FormatedVector& setSuffix(std::string suffix_);
        FormatedVector& setSeparator(std::string separator_);
    };

    template<Vector T>
    FormatedVector<T>::FormatedVector(const T& data_)
            : data(data_)
            , prefix("(")
            , suffix(")")
            , separator(", ") {}

    template<Vector U>
    std::ostream& operator<<(std::ostream& os, const FormatedVector<U>& v) {
        const U& data = v.data;
        os << v.prefix;
        size_t length = data.getLength();
        if (length > 0) {
            --length;
            for (size_t i = 0; i < length; ++i)
                os << data.calc(i) << v.separator;
            os << data.calc(length);
        }
        os << v.suffix;
        return os;
    }

    template<Vector T>
    FormatedVector<T>& FormatedVector<T>::toFormatMMA() {
        setPrefix("{");
        setSuffix("}");
        setSeparator(",");
        return *this;
    }

    template<Vector T>
    FormatedVector<T>& FormatedVector<T>::setPrefix(std::string prefix_) {
        prefix = std::move(prefix_);
        return *this;
    }

    template<Vector T>
    FormatedVector<T>& FormatedVector<T>::setSuffix(std::string suffix_) {
        suffix = std::move(suffix_);
        return *this;
    }

    template<Vector T>
    FormatedVector<T>& FormatedVector<T>::setSeparator(std::string separator_) {
        separator = std::move(separator_);
        return *this;
    }
}
