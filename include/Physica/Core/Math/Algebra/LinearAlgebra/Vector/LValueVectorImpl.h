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

namespace Physica::Core {
    template<class Derived>
    template<class PacketType>
    void LValueVector<Derived>::writePacket(size_t index, const PacketType packet) {
        using TrivialType = typename ScalarType::TrivialType;
        TrivialType buffer[PacketType::size()];
        packet.store(buffer);
        for (int i = 0; i < PacketType::size(); ++i, ++index)
            (*this)[index] = buffer[i];
    }

    template<class Derived>
    template<class PacketType>
    void LValueVector<Derived>::writePacketPartial(size_t index, const PacketType packet) {
        using TrivialType = typename ScalarType::TrivialType;
        TrivialType buffer[PacketType::size()];
        packet.store(buffer);
        for (int i = 0; index < Base::getLength(); ++i, ++index)
            (*this)[index] = buffer[i];
    }
}
