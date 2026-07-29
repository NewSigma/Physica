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

#include <H5Tpublic.h>
#include "H5ID.h"

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
    private:
        explicit H5Type(hid_t hid);
    };

    template<class T>
    auto H5Type::get() noexcept -> This {
        if constexpr (std::is_same_v<T, int8_t>)
            return This(H5T_NATIVE_INT8);
        else if constexpr (std::is_same_v<T, int16_t>)
            return This(H5T_NATIVE_INT16);
        else if constexpr (std::is_same_v<T, int32_t>)
            return This(H5T_NATIVE_INT32);
        else if constexpr (std::is_same_v<T, int64_t>)
            return This(H5T_NATIVE_INT64);
        else if constexpr (std::is_same_v<T, uint8_t>)
            return This(H5T_NATIVE_UINT8);
        else if constexpr (std::is_same_v<T, uint16_t>)
            return This(H5T_NATIVE_UINT16);
        else if constexpr (std::is_same_v<T, uint32_t>)
            return This(H5T_NATIVE_UINT32);
        else if constexpr (std::is_same_v<T, uint64_t>)
            return This(H5T_NATIVE_UINT64);
        else if constexpr (std::is_same_v<T, bool>)
            return This(H5T_NATIVE_HBOOL);
        else if constexpr (std::is_same_v<T, char>)
            return This(H5T_NATIVE_CHAR);
        else if constexpr (std::is_same_v<T, signed char>)
            return This(H5T_NATIVE_SCHAR);
        else if constexpr (std::is_same_v<T, unsigned char>)
            return This(H5T_NATIVE_UCHAR);
        else if constexpr (std::is_same_v<T, short>)
            return This(H5T_NATIVE_SHORT);
        else if constexpr (std::is_same_v<T, unsigned short>)
            return This(H5T_NATIVE_USHORT);
        else if constexpr (std::is_same_v<T, int>)
            return This(H5T_NATIVE_INT);
        else if constexpr (std::is_same_v<T, unsigned int>)
            return This(H5T_NATIVE_UINT);
        else if constexpr (std::is_same_v<T, long>)
            return This(H5T_NATIVE_LONG);
        else if constexpr (std::is_same_v<T, unsigned long>)
            return This(H5T_NATIVE_ULONG);
        else if constexpr (std::is_same_v<T, long long>)
            return This(H5T_NATIVE_LLONG);
        else if constexpr (std::is_same_v<T, unsigned long long>)
            return This(H5T_NATIVE_ULLONG);
        else if constexpr (std::is_same_v<T, float>)
            return This(H5T_NATIVE_FLOAT);
        else if constexpr (std::is_same_v<T, double>)
            return This(H5T_NATIVE_DOUBLE);
        else if constexpr (std::is_same_v<T, long double>)
            return This(H5T_NATIVE_LDOUBLE);
        else if constexpr (std::is_enum_v<T>)
            return get<std::underlying_type_t<T>>();
        else {
            using M = Traits<T>::MachineType;
            static_assert(!std::same_as<T, M>, "[Error]: Bad machine type");
            return get<M>();
        }
    }

    template<class T>
    auto H5Type::compound() -> This {
        return This(H5Tcreate(H5T_COMPOUND, sizeof(T)));
    }
}
