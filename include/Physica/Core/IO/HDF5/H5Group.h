/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    class PHYSICA_API H5Group : public H5::Group, public H5Loc {
        using Base = H5::Group;
        using Location = H5Loc;
    public:
        H5Group(H5::Group group) : Base(group) {}
        H5Group(const H5Group&) = default;
        H5Group(H5Group&&) noexcept = default;
        ~H5Group() = default;
        /* Operators */
        H5Group& operator=(H5Group& obj) = default;
        H5Group& operator=(H5Group&&) noexcept = delete;
        /* Operations */
        using Location::createDataSet;
        using Location::openDataSet;
        using Location::openGroup;

        template<class T>
        const H5::Attribute readAttr(const std::string& name, T& value) const;
        template<class T>
        H5::Attribute writeAttr(const std::string& name, T value);
    private:
        using Base::createDataSet;
        using Base::openDataSet;
        template<class T>
        [[nodiscard]] static const H5::PredType& getPredType();
    };

    template<class T>
    const H5::Attribute H5Group::readAttr(const std::string& name, T& value) const {
        constexpr bool IsArray = std::is_array<T>::value;
        constexpr size_t NumElem = IsArray ? std::extent<T>::value : 1;
        static_assert(!IsArray || std::rank<T>::value == 1, "[Error]: High dim array is not supported");
        static_assert(NumElem > 0, "[Error]: Bad array size");

        const auto type = getPredType<T>();
        const auto space = H5DataSpace<1>(NumElem);
        H5::Attribute attr;
        if (Base::attrExists(name.c_str()))
            attr = Base::openAttribute(name.c_str());
        else
            attr = Base::createAttribute(name.c_str(), type, space);
        attr.read(type, &value);
        return attr;
    }

    template<class T>
    H5::Attribute H5Group::writeAttr(const std::string& name, T value) {
        constexpr bool IsArray = std::is_array<T>::value;
        constexpr size_t NumElem = IsArray ? std::extent<T>::value : 1;
        static_assert(!IsArray || std::rank<T>::value == 1, "[Error]: High dim array is not supported");
        static_assert(NumElem > 0, "[Error]: Bad array size");

        const auto type = getPredType<T>();
        const auto space = H5DataSpace<1>(NumElem);
        H5::Attribute attr;
        if (Base::attrExists(name.c_str()))
            attr = Base::openAttribute(name.c_str());
        else
            attr = Base::createAttribute(name.c_str(), type, space);
        attr.write(type, &value);
        return attr;
    }

    template<class T>
    const H5::PredType& H5Group::getPredType() {
        if constexpr (std::is_same<T, int8_t>::value)
            return H5::PredType::NATIVE_INT8;
        else if constexpr (std::is_same<T, int16_t>::value)
            return H5::PredType::NATIVE_INT16;
        else if constexpr (std::is_same<T, int32_t>::value)
            return H5::PredType::NATIVE_INT32;
        else if constexpr (std::is_same<T, int64_t>::value)
            return H5::PredType::NATIVE_INT64;
        else if constexpr (std::is_same<T, uint8_t>::value)
            return H5::PredType::NATIVE_UINT8;
        else if constexpr (std::is_same<T, uint16_t>::value)
            return H5::PredType::NATIVE_UINT16;
        else if constexpr (std::is_same<T, uint32_t>::value)
            return H5::PredType::NATIVE_UINT32;
        else if constexpr (std::is_same<T, uint64_t>::value)
            return H5::PredType::NATIVE_UINT64;
        else if constexpr (std::is_same<T, float>::value)
            return H5::PredType::NATIVE_FLOAT;
        else if constexpr (std::is_same<T, Real<Float32>>::value)
            return H5::PredType::NATIVE_FLOAT;
        else if constexpr (std::is_same<T, double>::value)
            return H5::PredType::NATIVE_DOUBLE;
        else {
            static_assert(std::is_same<T, Real<Float64>>::value, "[Error]: Not implemented");
            return H5::PredType::NATIVE_DOUBLE;
        }
    }
}
