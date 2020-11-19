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
#ifndef PHYSICA_RINGBUFFER_H
#define PHYSICA_RINGBUFFER_H

#include <cstddef>

namespace Physica::Utils {
    class RingBuffer {
        /*!
         * Buffer of data bytes.
         */
        char* buffer;
        /*!
         * Size of buffer.
         */
        size_t size;
        /*!
         * Directs to the next position of buffer to be read.
         * It is always behind of bufferReader.
         */
        char* bufferReader;
        /*!
         * Directs to the next available position of buffer to write.
         * It is always in front of bufferReader.
         */
        char* bufferWriter;
    public:
        explicit RingBuffer(size_t size);
        RingBuffer(const RingBuffer& ring);
        RingBuffer(RingBuffer&& ring) noexcept;
        ~RingBuffer();
        /* Operators */
        RingBuffer& operator=(const RingBuffer& ring);
        RingBuffer& operator=(RingBuffer&& ring) noexcept;
        /* Operations */
        template<typename T>
        inline void write(T t);

        template<typename T>
        inline void read(T* t);
        /* Getters */
        [[nodiscard]] bool isEmpty() const noexcept { return bufferReader == bufferWriter; }
    private:
        void writeBytes(const char* src, size_t bytes);
        void readBytes(char* dest, size_t bytes);
    };
    /*!
     * Write T to the buffer.
     *
     * \tparam T
     * Arbitrary type.
     *
     * \param t
     * Data of t will be stored to the buffer.
     */
    template<typename T>
    inline void RingBuffer::write(T t) {
        writeBytes(reinterpret_cast<const char*>(&t), sizeof(T));
    }
    /*!
     * Read a T from the buffer and store it to t.
     *
     * \tparam T
     * Arbitrary type.
     *
     * \param t
     * The data of T will be save to the position directed by t.
     */
    template<typename T>
    inline void RingBuffer::read(T* t) {
        readBytes(reinterpret_cast<char*>(t), sizeof(T));
    }
}

#endif
