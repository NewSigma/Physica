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

namespace Physica {
    class PHYSICA_API MPIContext final {
        using This = MPIContext;
        enum class Dtype : int8_t {
            Int8,
            Int16,
            Int32,
            Int64,
            UInt8,
            UInt16,
            UInt32,
            UInt64,
            Bool,
            Char,
            SignedChar,
            UnsignedChar,
            Short,
            UnsignedShort,
            Int,
            UnsignedInt,
            Long,
            UnsignedLong,
            LongLong,
            UnsignedLongLong,
            Float,
            Double,
            LongDouble
        };
    public:
        using comm_handle = Handle<HandleType::MPI_Comm>;
        using dtype_handle = Handle<HandleType::MPI_Dtype>;

        static const comm_handle World;
    public:
        MPIContext(const This&) = delete;
        MPIContext(This&&) noexcept = delete;
        ~MPIContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static bool initialized() noexcept;
        [[nodiscard]] static int getNumProcess() noexcept;
        [[nodiscard]] static int getProcessID() noexcept;

        static void send(int from, int to, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void sendrecv(int send_to, int recv_from, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void bcast(int root, void* data, int count, dtype_handle dtype, comm_handle comm = World);
        static void wait(comm_handle comm = World);

        template<class T>
        [[nodiscard]] static dtype_handle dtype() noexcept;
    private:
        MPIContext() noexcept;
        /* Static members */
        [[nodiscard]] static dtype_handle dtype_scalar(Dtype type) noexcept;
    };

    template<class T>
    auto MPIContext::dtype() noexcept -> dtype_handle {
        if constexpr (std::is_same_v<T, int8_t>)
            return dtype_scalar(Dtype::Int8);
        else if constexpr (std::is_same_v<T, int16_t>)
            return dtype_scalar(Dtype::Int16);
        else if constexpr (std::is_same_v<T, int32_t>)
            return dtype_scalar(Dtype::Int32);
        else if constexpr (std::is_same_v<T, int64_t>)
            return dtype_scalar(Dtype::Int64);
        else if constexpr (std::is_same_v<T, uint8_t>)
            return dtype_scalar(Dtype::UInt8);
        else if constexpr (std::is_same_v<T, uint16_t>)
            return dtype_scalar(Dtype::UInt16);
        else if constexpr (std::is_same_v<T, uint32_t>)
            return dtype_scalar(Dtype::UInt32);
        else if constexpr (std::is_same_v<T, uint64_t>)
            return dtype_scalar(Dtype::UInt64);
        else if constexpr (std::is_same_v<T, bool>)
            return dtype_scalar(Dtype::Bool);
        else if constexpr (std::is_same_v<T, char>)
            return dtype_scalar(Dtype::Char);
        else if constexpr (std::is_same_v<T, signed char>)
            return dtype_scalar(Dtype::SignedChar);
        else if constexpr (std::is_same_v<T, unsigned char>)
            return dtype_scalar(Dtype::UnsignedChar);
        else if constexpr (std::is_same_v<T, short>)
            return dtype_scalar(Dtype::Short);
        else if constexpr (std::is_same_v<T, unsigned short>)
            return dtype_scalar(Dtype::UnsignedShort);
        else if constexpr (std::is_same_v<T, int>)
            return dtype_scalar(Dtype::Int);
        else if constexpr (std::is_same_v<T, unsigned int>)
            return dtype_scalar(Dtype::UnsignedInt);
        else if constexpr (std::is_same_v<T, long>)
            return dtype_scalar(Dtype::Long);
        else if constexpr (std::is_same_v<T, unsigned long>)
            return dtype_scalar(Dtype::UnsignedLong);
        else if constexpr (std::is_same_v<T, long long>)
            return dtype_scalar(Dtype::LongLong);
        else if constexpr (std::is_same_v<T, unsigned long long>)
            return dtype_scalar(Dtype::UnsignedLongLong);
        else if constexpr (std::is_same_v<T, float>)
            return dtype_scalar(Dtype::Float);
        else if constexpr (std::is_same_v<T, double>)
            return dtype_scalar(Dtype::Double);
        else if constexpr (std::is_same_v<T, long double>)
            return dtype_scalar(Dtype::LongDouble);
        else if constexpr (std::is_enum_v<T>)
            return dtype<std::underlying_type_t<T>>();
        else
            return T::dtype_mpi();
    }
}
