/*
 * Copyright 2020 WeiBo He.
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
#include "Physica/Utils/RingBuffer.h"

namespace Physica::Utils {
    RingBuffer::RingBuffer(size_t size)
            : buffer(new char[size])
            , size(size) {
        bufferReader = bufferWriter = buffer;
    }

    RingBuffer::~RingBuffer() {
        delete[] buffer;
    }

    RingBuffer::RingBuffer(const RingBuffer& ring)
            : buffer(new char[ring.size])
            , size(ring.size)
            , bufferReader(ring.bufferReader)
            , bufferWriter(ring.bufferWriter) {}

    RingBuffer::RingBuffer(RingBuffer&& ring) noexcept
            : buffer(ring.buffer)
            , size(ring.size)
            , bufferReader(ring.bufferReader)
            , bufferWriter(ring.bufferWriter) {
        ring.buffer = nullptr;
    }

    RingBuffer& RingBuffer::operator=(const RingBuffer& ring) {
        if(this != &ring) {
            this->~RingBuffer();
            size = ring.size;
            buffer = new char[size];
            bufferReader = ring.bufferReader;
            bufferWriter = ring.bufferWriter;
        }
        return *this;
    }

    RingBuffer& RingBuffer::operator=(RingBuffer&& ring) noexcept {
        this->~RingBuffer();
        buffer = ring.buffer;
        ring.buffer = nullptr;
        size = ring.size;
        bufferReader = ring.bufferReader;
        bufferWriter = ring.bufferWriter;
        return *this;
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
        assert(bytes < size);
        size_t leftSpace = size - (bufferWriter - buffer);
        if(bytes < leftSpace) {
            while(bufferReader > bufferWriter && (bufferReader - bufferWriter) < bytes)
                std::this_thread::yield();
            memcpy(bufferWriter, src, bytes);
            bufferWriter = bytes == leftSpace ? buffer : (bufferWriter + bytes);
        }
        else {
            size_t leftBytes = bytes - leftSpace;
            while(bufferReader > bufferWriter || (bufferReader - buffer) < leftBytes)
                std::this_thread::yield();
            memcpy(bufferWriter, src, leftSpace);
            memcpy(buffer, src + leftSpace, leftBytes);
            bufferWriter = buffer + leftBytes;
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
        size_t leftSpace = size - (bufferReader - buffer);
        if(bytes <= leftSpace) {
            while(bufferReader == bufferWriter)
                std::this_thread::yield();
            memcpy(dest, bufferReader, bytes);
            bufferReader = bytes == leftSpace ? buffer : (bufferReader + bytes);
        }
        else {
            size_t leftBytes = bytes - leftSpace;
            while(bufferReader == bufferWriter)
                std::this_thread::yield();
            memcpy(dest, bufferReader, leftSpace);
            memcpy(dest + leftSpace, buffer, leftBytes);
            bufferReader = buffer + leftBytes;
        }
    }
}
