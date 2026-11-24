/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Array.h"

namespace Physica {
    class PHYSICA_API RingBuffer {
        using This = RingBuffer;

        Array<char> buffer;
        char* bufferReader;
        char* bufferWriter;
    public:
        explicit RingBuffer(size_t size);
        RingBuffer(const This&) = default;
        RingBuffer(This&&) noexcept = default;
        ~RingBuffer() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<typename T> void write(T t) noexcept;
        template<typename T> void read(T* __restrict t) noexcept;
        template<typename T> void cread(T* __restrict t, size_t bias) const noexcept;
        void writeBytes(const char* __restrict src, size_t bytes) noexcept;
        void readBytes(char* __restrict dest, size_t bytes) noexcept;
        void creadBytes(char* __restrict dest, size_t bytes, size_t bias) const noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool empty() const noexcept { return bufferReader == bufferWriter; }
        [[nodiscard]] size_t getSize() const noexcept { return buffer.getLength(); }
    };

    template<typename T>
    void RingBuffer::write(T t) noexcept {
        static_assert(std::is_trivially_copyable<T>::value);
        writeBytes(reinterpret_cast<const char*>(&t), sizeof(T));
    }

    template<typename T>
    void RingBuffer::read(T* __restrict t) noexcept {
        static_assert(std::is_trivially_copyable<T>::value);
        readBytes(reinterpret_cast<char*>(t), sizeof(T));
    }

    template<typename T>
    void RingBuffer::cread(T* __restrict t, size_t bias) const noexcept {
        static_assert(std::is_trivially_copyable<T>::value);
        creadBytes(reinterpret_cast<char*>(t), sizeof(T), bias);
    }
}
