/*
 * Copyright 2026 Weibo He.
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

#include <cstdint>
#include <type_traits>
#include <H5Tpublic.h>
#include "Physica/Core/Utils/PrimitiveType.h"
#include "Mixin/H5ID.h"

namespace Physica {
    class PHYSICA_API H5Type : public H5ID {
        using This = H5Type;
    public:
        H5Type(const This&) noexcept = default;
        H5Type(This&&) noexcept = default;
        ~H5Type() noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        void insert(const char* name, size_t offset, const H5Type& memberType) noexcept;
        /* Getters */
        [[nodiscard]] bool isCompound() const noexcept;
        /* Static members */
        template<class T>
        [[nodiscard]] static This get() noexcept;
        template<class T>
        [[nodiscard]] static This compound();
        [[nodiscard]] constexpr static IdentifierType itype() noexcept { return IdentifierType::Datatype; }
    private:
        explicit H5Type(H5ID id) noexcept;
        explicit H5Type(hid_t hid) noexcept;
        /* Static members */
        [[nodiscard]] static This getPrimitiveType(PrimitiveType type) noexcept;

        friend class H5ID;
    };

    template<class T>
    auto H5Type::get() noexcept -> This {
        if constexpr (std::is_same_v<T, int8_t>)
            return getPrimitiveType(PrimitiveType::Int8);
        else if constexpr (std::is_same_v<T, int16_t>)
            return getPrimitiveType(PrimitiveType::Int16);
        else if constexpr (std::is_same_v<T, int32_t>)
            return getPrimitiveType(PrimitiveType::Int32);
        else if constexpr (std::is_same_v<T, int64_t>)
            return getPrimitiveType(PrimitiveType::Int64);
        else if constexpr (std::is_same_v<T, uint8_t>)
            return getPrimitiveType(PrimitiveType::UInt8);
        else if constexpr (std::is_same_v<T, uint16_t>)
            return getPrimitiveType(PrimitiveType::UInt16);
        else if constexpr (std::is_same_v<T, uint32_t>)
            return getPrimitiveType(PrimitiveType::UInt32);
        else if constexpr (std::is_same_v<T, uint64_t>)
            return getPrimitiveType(PrimitiveType::UInt64);
        else if constexpr (std::is_same_v<T, bool>)
            return getPrimitiveType(PrimitiveType::Bool);
        else if constexpr (std::is_same_v<T, char>)
            return getPrimitiveType(PrimitiveType::Char);
        else if constexpr (std::is_same_v<T, signed char>)
            return getPrimitiveType(PrimitiveType::SignedChar);
        else if constexpr (std::is_same_v<T, unsigned char>)
            return getPrimitiveType(PrimitiveType::UnsignedChar);
        else if constexpr (std::is_same_v<T, short>)
            return getPrimitiveType(PrimitiveType::Short);
        else if constexpr (std::is_same_v<T, unsigned short>)
            return getPrimitiveType(PrimitiveType::UnsignedShort);
        else if constexpr (std::is_same_v<T, int>)
            return getPrimitiveType(PrimitiveType::Int);
        else if constexpr (std::is_same_v<T, unsigned int>)
            return getPrimitiveType(PrimitiveType::UnsignedInt);
        else if constexpr (std::is_same_v<T, long>)
            return getPrimitiveType(PrimitiveType::Long);
        else if constexpr (std::is_same_v<T, unsigned long>)
            return getPrimitiveType(PrimitiveType::UnsignedLong);
        else if constexpr (std::is_same_v<T, long long>)
            return getPrimitiveType(PrimitiveType::LongLong);
        else if constexpr (std::is_same_v<T, unsigned long long>)
            return getPrimitiveType(PrimitiveType::UnsignedLongLong);
        else if constexpr (std::is_same_v<T, float>)
            return getPrimitiveType(PrimitiveType::Float);
        else if constexpr (std::is_same_v<T, double>)
            return getPrimitiveType(PrimitiveType::Double);
        else if constexpr (std::is_same_v<T, long double>)
            return getPrimitiveType(PrimitiveType::LongDouble);
        else if constexpr (std::is_enum_v<T>)
            return get<std::underlying_type_t<T>>();
        else {
            using M = Traits<T>::MachineType;
            static_assert(!std::is_same_v<T, M>, "[Error]: Bad machine type");
            return get<M>();
        }
    }

    template<class T>
    auto H5Type::compound() -> This {
        return This(H5Tcreate(H5T_COMPOUND, sizeof(T)));
    }
}
