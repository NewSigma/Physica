/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstddef>
#include "Physica/Macro.h"

namespace Physica {
    /**
     * BitArray is similar to bool[] but it has a better use of space.
     */
    class PHYSICA_API BitArray {
        using This = BitArray;
        //Every bit of arr stands for a bool.
        unsigned char* arr = nullptr;
        size_t bitCount = 0;
    public:
        BitArray() = default;
        explicit BitArray(size_t bitCount_);
        BitArray(size_t bitCount_, bool initial);
        BitArray(const BitArray& array);
        BitArray(BitArray&& array) noexcept;
        ~BitArray();
        /* Operators */
        BitArray& operator=(const BitArray& array);
        BitArray& operator=(BitArray&& array) noexcept;
        bool operator[](size_t s) const;
        BitArray operator&(const BitArray& array) const;
        BitArray operator|(const BitArray& array) const;
        BitArray operator~() const;
        /* Operations */
        void setBit(size_t s, bool b);
        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept;
    private:
        BitArray(unsigned char* arr, size_t bitCount);
    };

    constexpr size_t BitArray::getLength() const noexcept {
        //((bitCount >> 3U) + (bitCount & 8U)) is the upper approximation of (bitCount / 8).
        return (bitCount >> 3U) + (bitCount & 8U);
    }
}
