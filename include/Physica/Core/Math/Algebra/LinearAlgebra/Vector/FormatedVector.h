/*
 * Copyright 2022 WeiBo He.
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

#include "RValueVector.h"

namespace Physica::Core {
    template<class VectorType>
    class FormatedVector {
        const RValueVector<VectorType>& data;
        std::string prefix;
        std::string suffix;
        std::string separator;
    public:
        FormatedVector(const RValueVector<VectorType>& data_);
        /* Operators */
        template<class T>
        friend std::ostream& operator<<(std::ostream& os, const FormatedVector<T>& v);
        /* Setters */
        FormatedVector& setPrefix(std::string prefix_);
        FormatedVector& setSuffix(std::string suffix_);
        FormatedVector& setSeparator(std::string separator_);
    };

    template<class VectorType>
    FormatedVector<VectorType>::FormatedVector(const RValueVector<VectorType>& data_)
            : data(data_)
            , prefix("(")
            , suffix(")")
            , separator(", ") {}

    template<class VectorType>
    std::ostream& operator<<(std::ostream& os, const FormatedVector<VectorType>& v) {
        const RValueVector<VectorType>& data = v.data;
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

    template<class VectorType>
    FormatedVector<VectorType>& FormatedVector<VectorType>::setPrefix(std::string prefix_) {
        prefix = std::move(prefix_);
        return *this;
    }

    template<class VectorType>
    FormatedVector<VectorType>& FormatedVector<VectorType>::setSuffix(std::string suffix_) {
        suffix = std::move(suffix_);
        return *this;
    }

    template<class VectorType>
    FormatedVector<VectorType>& FormatedVector<VectorType>::setSeparator(std::string separator_) {
        separator = std::move(separator_);
        return *this;
    }
}
