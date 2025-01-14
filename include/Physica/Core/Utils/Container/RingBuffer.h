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
#pragma once

#include "Array.h"

namespace Physica::Core {
    class PHYSICA_API RingBuffer {
        using This = RingBuffer;

        Array<char> buffer;
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
        RingBuffer(const This&) = default;
        RingBuffer(This&&) noexcept = default;
        ~RingBuffer() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<typename T> inline void write(T t);
        template<typename T> inline void read(T* t);
        template<typename T> inline void cread(T* t, size_t bias) const;
        void writeBytes(const char* src, size_t bytes);
        void readBytes(char* dest, size_t bytes);
        void creadBytes(char* dest, size_t bytes, size_t bias) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool isEmpty() const noexcept { return bufferReader == bufferWriter; }
        [[nodiscard]] size_t getSize() const noexcept { return buffer.getLength(); }
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

    template<typename T>
    inline void RingBuffer::cread(T* t, size_t bias) const {
        creadBytes(reinterpret_cast<char*>(t), sizeof(T), bias);
    }
}
