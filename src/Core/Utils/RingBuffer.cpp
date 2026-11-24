/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cassert>
#include <cstring>
#include <thread>
#include "Physica/Core/Utils/Container/RingBuffer.h"

using namespace Physica;

RingBuffer::RingBuffer(size_t size) : buffer(size) {
    bufferReader = bufferWriter = buffer.data();
}

void RingBuffer::writeBytes(const char*  __restrict src, size_t bytes) noexcept {
    assert(bytes < getSize());
    const size_t leftSpace = getSize() - (bufferWriter - buffer.data());
    if (bytes < leftSpace) {
        while (bufferReader > bufferWriter && static_cast<size_t>(bufferReader - bufferWriter) < bytes)
            std::this_thread::yield();
        memcpy(bufferWriter, src, bytes);
        bufferWriter = bytes == leftSpace ? buffer.data() : (bufferWriter + bytes);
    }
    else {
        size_t leftBytes = bytes - leftSpace;
        while (bufferReader > bufferWriter || static_cast<size_t>(bufferReader - buffer.data()) < leftBytes)
            std::this_thread::yield();
        memcpy(bufferWriter, src, leftSpace);
        memcpy(buffer.data(), src + leftSpace, leftBytes);
        bufferWriter = buffer.data() + leftBytes;
    }
}

void RingBuffer::readBytes(char*  __restrict dest, size_t bytes) noexcept {
    const size_t leftSpace = getSize() - (bufferReader - buffer.data());
    if (bytes <= leftSpace) {
        while (bufferReader == bufferWriter)
            std::this_thread::yield();
        memcpy(dest, bufferReader, bytes);
        bufferReader = bytes == leftSpace ? buffer.data() : (bufferReader + bytes);
    }
    else {
        size_t leftBytes = bytes - leftSpace;
        while (bufferReader == bufferWriter)
            std::this_thread::yield();
        memcpy(dest, bufferReader, leftSpace);
        memcpy(dest + leftSpace, buffer.data(), leftBytes);
        bufferReader = buffer.data() + leftBytes;
    }
}

void RingBuffer::creadBytes(char*  __restrict dest, size_t bytes, size_t bias) const noexcept {
    const char* startPos = buffer.data() + (bufferReader - buffer.data() + bias) % getSize();
    const size_t leftSpace = getSize() - (startPos - buffer.data());
    if (bytes <= leftSpace) {
        while (startPos <= bufferWriter && startPos + bytes > bufferWriter)
            std::this_thread::yield();
        memcpy(dest, startPos, bytes);
    }
    else {
        size_t leftBytes = bytes - leftSpace;
        while (buffer.data() + leftBytes > bufferWriter)
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
