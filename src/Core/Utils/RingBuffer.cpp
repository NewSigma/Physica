/*
 * Copyright 2020-2025 Weibo He.
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
#include <cassert>
#include <cstring>
#include <thread>
#include "Physica/Core/Utils/Container/RingBuffer.h"

namespace Physica {
    RingBuffer::RingBuffer(size_t size) : buffer(size) {
        bufferReader = bufferWriter = buffer.data();
    }
    /*!
     * Read data from src and write bytes bytes to the buffer.
     *
     * \param src
     * Read bytes from it.
     *
     * \param bytes
     * Number of bytes to be read and write.
     */
    void RingBuffer::writeBytes(const char* src, size_t bytes) {
        assert(bytes < getSize());
        const size_t leftSpace = getSize() - (bufferWriter - buffer.data());
        if(bytes < leftSpace) {
            while(bufferReader > bufferWriter && static_cast<size_t>(bufferReader - bufferWriter) < bytes)
                std::this_thread::yield();
            memcpy(bufferWriter, src, bytes);
            bufferWriter = bytes == leftSpace ? buffer.data() : (bufferWriter + bytes);
        }
        else {
            size_t leftBytes = bytes - leftSpace;
            while(bufferReader > bufferWriter || static_cast<size_t>(bufferReader - buffer.data()) < leftBytes)
                std::this_thread::yield();
            memcpy(bufferWriter, src, leftSpace);
            memcpy(buffer.data(), src + leftSpace, leftBytes);
            bufferWriter = buffer.data() + leftBytes;
        }
    }
    /*!
     * Read data from buffer and write bytes bytes to the dest.
     *
     * \param dest
     * Save bytes to it.
     *
     * \param bytes
     * Number of bytes to be read and write.
     */
    void RingBuffer::readBytes(char* dest, size_t bytes) {
        const size_t leftSpace = getSize() - (bufferReader - buffer.data());
        if(bytes <= leftSpace) {
            while(bufferReader == bufferWriter)
                std::this_thread::yield();
            memcpy(dest, bufferReader, bytes);
            bufferReader = bytes == leftSpace ? buffer.data() : (bufferReader + bytes);
        }
        else {
            size_t leftBytes = bytes - leftSpace;
            while(bufferReader == bufferWriter)
                std::this_thread::yield();
            memcpy(dest, bufferReader, leftSpace);
            memcpy(dest + leftSpace, buffer.data(), leftBytes);
            bufferReader = buffer.data() + leftBytes;
        }
    }

    void RingBuffer::creadBytes(char* dest, size_t bytes, size_t bias) const {
        const char* startPos = buffer.data() + (bufferReader - buffer.data() + bias) % getSize();
        const size_t leftSpace = getSize() - (startPos - buffer.data());
        if(bytes <= leftSpace) {
            while(startPos <= bufferWriter && startPos + bytes > bufferWriter)
                std::this_thread::yield();
            memcpy(dest, startPos, bytes);
        }
        else {
            size_t leftBytes = bytes - leftSpace;
            while(buffer.data() + leftBytes > bufferWriter)
                std::this_thread::yield();
            memcpy(dest, startPos, leftSpace);
            memcpy(dest + leftSpace, buffer.data(), leftBytes);
        }
    }

    void RingBuffer::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        buffer.swap(obj.buffer);
        std::swap(bufferReader, obj.bufferReader);
        std::swap(bufferWriter, obj.bufferWriter);
    }
}
