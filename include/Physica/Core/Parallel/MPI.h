/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Handle.h"
#include "Physica/Core/Utils/PrimitiveType.h"

namespace Physica {
    class PHYSICA_API MPI final {
    public:
        using comm_handle = Handle<HandleType::MPI_Comm>;
        using dtype_handle = Handle<HandleType::MPI_Dtype>;

        static const comm_handle World;
    public:
        MPI(const MPI&) = delete;
        MPI(MPI&&) noexcept = delete;
        ~MPI();
        /* Operators */
        MPI& operator=(const MPI&) = delete;
        MPI& operator=(MPI&&) noexcept = delete;
        /* Static members */
        [[nodiscard]] static MPI& getInstance() noexcept;
        [[nodiscard]] static bool initialized() noexcept;
        [[nodiscard]] static int getNumProcess() noexcept;
        [[nodiscard]] static int getProcessID() noexcept;

        static void send(int to, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void recv(int from, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void pass(int from, int to, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void sendrecv(int send_to, int recv_from, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void bcast(int root, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void wait(comm_handle comm = World);

        template<class T>
        [[nodiscard]] static dtype_handle dtype() noexcept;
    private:
        MPI() noexcept;
        /* Static memebers */
        [[nodiscard]] static dtype_handle dtype_primitive(PrimitiveType type) noexcept;
    };

    template<class T>
    auto MPI::dtype() noexcept -> dtype_handle {
        if constexpr (std::is_same_v<T, int8_t>)
            return dtype_primitive(PrimitiveType::Int8);
        else if constexpr (std::is_same_v<T, int16_t>)
            return dtype_primitive(PrimitiveType::Int16);
        else if constexpr (std::is_same_v<T, int32_t>)
            return dtype_primitive(PrimitiveType::Int32);
        else if constexpr (std::is_same_v<T, int64_t>)
            return dtype_primitive(PrimitiveType::Int64);
        else if constexpr (std::is_same_v<T, uint8_t>)
            return dtype_primitive(PrimitiveType::UInt8);
        else if constexpr (std::is_same_v<T, uint16_t>)
            return dtype_primitive(PrimitiveType::UInt16);
        else if constexpr (std::is_same_v<T, uint32_t>)
            return dtype_primitive(PrimitiveType::UInt32);
        else if constexpr (std::is_same_v<T, uint64_t>)
            return dtype_primitive(PrimitiveType::UInt64);
        else if constexpr (std::is_same_v<T, bool>)
            return dtype_primitive(PrimitiveType::Bool);
        else if constexpr (std::is_same_v<T, char>)
            return dtype_primitive(PrimitiveType::Char);
        else if constexpr (std::is_same_v<T, signed char>)
            return dtype_primitive(PrimitiveType::SignedChar);
        else if constexpr (std::is_same_v<T, unsigned char>)
            return dtype_primitive(PrimitiveType::UnsignedChar);
        else if constexpr (std::is_same_v<T, short>)
            return dtype_primitive(PrimitiveType::Short);
        else if constexpr (std::is_same_v<T, unsigned short>)
            return dtype_primitive(PrimitiveType::UnsignedShort);
        else if constexpr (std::is_same_v<T, int>)
            return dtype_primitive(PrimitiveType::Int);
        else if constexpr (std::is_same_v<T, unsigned int>)
            return dtype_primitive(PrimitiveType::UnsignedInt);
        else if constexpr (std::is_same_v<T, long>)
            return dtype_primitive(PrimitiveType::Long);
        else if constexpr (std::is_same_v<T, unsigned long>)
            return dtype_primitive(PrimitiveType::UnsignedLong);
        else if constexpr (std::is_same_v<T, long long>)
            return dtype_primitive(PrimitiveType::LongLong);
        else if constexpr (std::is_same_v<T, unsigned long long>)
            return dtype_primitive(PrimitiveType::UnsignedLongLong);
        else if constexpr (std::is_same_v<T, float>)
            return dtype_primitive(PrimitiveType::Float);
        else if constexpr (std::is_same_v<T, double>)
            return dtype_primitive(PrimitiveType::Double);
        else if constexpr (std::is_same_v<T, long double>)
            return dtype_primitive(PrimitiveType::LongDouble);
        else if constexpr (std::is_enum_v<T>)
            return dtype<std::underlying_type_t<T>>();
        else
            return T::dtype_mpi();
    }
}
