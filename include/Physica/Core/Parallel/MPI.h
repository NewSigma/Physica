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
#include "Physica/Core/Utils/PrimitiveType.h"
#include "MPI/Request.h"

namespace Physica {
    class PHYSICA_API MPI final {
    public:
        using comm_handle = Handle<HandleType::MPI_Comm>;
        using dtype_handle = Handle<HandleType::MPI_Dtype>;
        using op_handle = Handle<HandleType::MPI_Op>;
        using request_handle = Request::Handle;

        enum class ReduceOp : int8_t {
            Sum,
            Prod,
            Max,
            Min,
            Land,
            Lor,
            BAnd,
            BOr,
            LXor,
            BXor,
        };

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
        [[nodiscard]] static int getNumRank() noexcept;
        [[nodiscard]] static int getRank() noexcept;

        static auto send(int to, const std::ranges::contiguous_range auto& buffer, comm_handle comm = World) -> Request;
        static auto recv(int from, std::ranges::contiguous_range auto& buffer, comm_handle comm = World) -> Request;
        static auto pass(int from, int to, std::ranges::contiguous_range auto& buffer, comm_handle comm = World) -> Request;
        static void sendrecv(int send_to, int recv_from, std::ranges::contiguous_range auto& buffer, comm_handle comm = World);
        static auto bcast(int root, std::ranges::contiguous_range auto& buffer, comm_handle comm = World) -> Request;
        static auto reduce(int to,
                           const std::ranges::contiguous_range auto& sendbuf,
                           std::ranges::contiguous_range auto& recvbuf,
                           ReduceOp op,
                           comm_handle comm = World) -> Request;
        static auto allreduce(const std::ranges::contiguous_range auto& sendbuf,
                              std::ranges::contiguous_range auto& recvbuf,
                              ReduceOp op,
                              comm_handle comm = World) -> Request;
        static auto gather(int to,
                           const std::ranges::contiguous_range auto& sendbuf,
                           std::ranges::contiguous_range auto& recvbuf,
                           comm_handle comm = World) -> Request;
        static auto scatter(int from,
                            const std::ranges::contiguous_range auto& sendbuf,
                            std::ranges::contiguous_range auto& recvbuf,
                            comm_handle comm = World) -> Request;
        static auto allgather(const std::ranges::contiguous_range auto& sendbuf,
                              std::ranges::contiguous_range auto& recvbuf,
                              comm_handle comm = World) -> Request;
        static void wait(comm_handle comm = World);

        template<class T>
        [[nodiscard]] static dtype_handle dtype() noexcept;
    private:
        MPI() noexcept;
        /* Static memebers */
        [[nodiscard]] static dtype_handle dtype_primitive(PrimitiveType type) noexcept;

        static auto send(int to, const void* data, int count, dtype_handle dtype, comm_handle comm) -> Request;
        static auto recv(int from, void* data, int count, dtype_handle dtype, comm_handle comm) -> Request;
        static auto pass(int from, int to, void* data, int count, dtype_handle dtype, comm_handle comm) -> Request;
        static void sendrecv(int send_to, int recv_from, void* data, int count, dtype_handle dtype, comm_handle comm);
        static auto bcast(int root, void* data, int count, dtype_handle dtype, comm_handle comm) -> Request;
        static auto reduce(int to, const void* sendbuf, void* recvbuf, int count, dtype_handle dtype, ReduceOp op, comm_handle comm) -> Request;
        static auto allreduce(const void* sendbuf, void* recvbuf, int count, dtype_handle dtype, ReduceOp op, comm_handle comm) -> Request;
        static auto gather(int to, const void* sendbuf, void* recvbuf, int count, dtype_handle dtype, comm_handle comm) -> Request;
        static auto scatter(int from, const void* sendbuf, void* recvbuf, int count, dtype_handle dtype, comm_handle comm) -> Request;
        static auto allgather(const void* sendbuf, void* recvbuf, int count, dtype_handle dtype, comm_handle comm) -> Request;
    };

    auto MPI::send(int to, const std::ranges::contiguous_range auto& buffer, comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(buffer)>;
        return send(to, std::ranges::data(buffer), std::ranges::size(buffer), dtype<T>(), comm);
    }

    auto MPI::recv(int from, std::ranges::contiguous_range auto& buffer, comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(buffer)>;
        return recv(from, std::ranges::data(buffer), std::ranges::size(buffer), dtype<T>(), comm);
    }

    auto MPI::pass(int from, int to, std::ranges::contiguous_range auto& buffer, comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(buffer)>;
        return pass(from, to, std::ranges::data(buffer), std::ranges::size(buffer), dtype<T>(), comm);
    }

    void MPI::sendrecv(int send_to, int recv_from, std::ranges::contiguous_range auto& buffer, comm_handle comm) {
        using T = std::ranges::range_value_t<decltype(buffer)>;
        sendrecv(send_to, recv_from, std::ranges::data(buffer), std::ranges::size(buffer), dtype<T>(), comm);
    }

    auto MPI::bcast(int root, std::ranges::contiguous_range auto& buffer, comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(buffer)>;
        return bcast(root, std::ranges::data(buffer), std::ranges::size(buffer), dtype<T>(), comm);
    }

    auto MPI::reduce(int to,
                     const std::ranges::contiguous_range auto& sendbuf,
                     std::ranges::contiguous_range auto& recvbuf,
                     ReduceOp op,
                     comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(sendbuf)>;
        using U = std::ranges::range_value_t<decltype(recvbuf)>;
        static_assert(std::same_as<T, U>, "[Error]: sendbuf and recvbuf must have the same value type");
        return reduce(to, std::ranges::data(sendbuf), std::ranges::data(recvbuf), std::ranges::size(sendbuf), dtype<T>(), op, comm);
    }

    auto MPI::allreduce(const std::ranges::contiguous_range auto& sendbuf,
                        std::ranges::contiguous_range auto& recvbuf,
                        ReduceOp op,
                        comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(sendbuf)>;
        using U = std::ranges::range_value_t<decltype(recvbuf)>;
        static_assert(std::same_as<T, U>, "[Error]: sendbuf and recvbuf must have the same value type");
        return allreduce(std::ranges::data(sendbuf), std::ranges::data(recvbuf), std::ranges::size(sendbuf), dtype<T>(), op, comm);
    }

    auto MPI::gather(int to,
                     const std::ranges::contiguous_range auto& sendbuf,
                     std::ranges::contiguous_range auto& recvbuf,
                     comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(sendbuf)>;
        using U = std::ranges::range_value_t<decltype(recvbuf)>;
        static_assert(std::same_as<T, U>, "[Error]: sendbuf and recvbuf must have the same value type");
        return gather(to, std::ranges::data(sendbuf), std::ranges::data(recvbuf), std::ranges::size(sendbuf), dtype<T>(), comm);
    }

    auto MPI::scatter(int from,
                      const std::ranges::contiguous_range auto& sendbuf,
                      std::ranges::contiguous_range auto& recvbuf,
                      comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(sendbuf)>;
        using U = std::ranges::range_value_t<decltype(recvbuf)>;
        static_assert(std::same_as<T, U>, "[Error]: sendbuf and recvbuf must have the same value type");
        return scatter(from, std::ranges::data(sendbuf), std::ranges::data(recvbuf), std::ranges::size(sendbuf), dtype<T>(), comm);
    }

    auto MPI::allgather(const std::ranges::contiguous_range auto& sendbuf,
                        std::ranges::contiguous_range auto& recvbuf,
                        comm_handle comm) -> Request {
        using T = std::ranges::range_value_t<decltype(sendbuf)>;
        using U = std::ranges::range_value_t<decltype(recvbuf)>;
        static_assert(std::same_as<T, U>, "[Error]: sendbuf and recvbuf must have the same value type");
        return allgather(std::ranges::data(sendbuf), std::ranges::data(recvbuf), std::ranges::size(sendbuf) / getNumRank(), dtype<T>(), comm);
    }

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
