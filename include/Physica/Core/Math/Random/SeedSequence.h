/*
 * Copyright 2025 Weibo He.
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

#include <algorithm>
#include <array>
#include <cstdint>

namespace Physica {
    /**
     * Based on randutils [1], refer to [2] for details.
     * We make it stateful for several generations.
     *
     * Reference:
     * [1] Melissa E. O'Neill, Random-Number Utilities (2015-2022); https://gist.github.com/imneme/540829265469e673d045
     * [2] Developing a seed_seq Alternative; http://www.pcg-random.org/posts/developing-a-seed_seq-alternative.html
     */
    template<size_t Count = 4, typename IntRep = uint32_t, size_t MixRound = 1 + (Count <= 2)>
    class SeedSequence {
        using This = SeedSequence<Count, IntRep, MixRound>;

        constexpr static uint32_t InitA = 0x43B0D7E5;
        constexpr static uint32_t MultA = 0x931E8875;
        constexpr static uint32_t InitB = 0x8B51F9DD;
        constexpr static uint32_t MultB = 0x58F38DED;
        constexpr static uint32_t MixMultL = 0xCA01F9DD;
        constexpr static uint32_t MixMultR = 0x4973F715;
        constexpr static uint32_t XShift = sizeof(IntRep) * 8 / 2;
    public:
        using result_type = IntRep;
    private:
        std::array<IntRep, Count> inits;
        std::array<IntRep, Count> mixer;
        uint32_t hasherA = InitA;
        uint32_t hasherB = InitB;
        int counter = 0;
    public:
        SeedSequence() = default;
        template<typename InputIter>
        SeedSequence(InputIter begin, InputIter end) noexcept;
        SeedSequence(std::initializer_list<IntRep> init) noexcept;
        SeedSequence(const This&) = default;
        SeedSequence(This&&) noexcept = default;
        ~SeedSequence() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<typename InputIter>
        void seed(InputIter begin, InputIter end) noexcept;

        template<std::integral T = IntRep>
        T generate() noexcept;
        template<typename RandomAccessIterator>
        void generate(RandomAccessIterator begin, RandomAccessIterator end) noexcept;

        template<typename OutputIterator>
        void param(OutputIterator dest) const noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t size() noexcept { return Count; }
        [[nodiscard]] const auto& getInits() const noexcept { return inits; }
    private:
        template<typename InputIter>
        void mix_entropy(InputIter begin, InputIter end) noexcept;
    };

    template<size_t Count, typename IntRep, size_t MixRound>
    template<typename InputIter>
    SeedSequence<Count, IntRep, MixRound>::SeedSequence(InputIter begin, InputIter end) noexcept {
        seed(begin, end);
    }

    template<size_t Count, typename IntRep, size_t MixRound>
    SeedSequence<Count, IntRep, MixRound>::SeedSequence(std::initializer_list<IntRep> init) noexcept : This(init.begin(), init.end()) {}

    template<size_t Count, typename IntRep, size_t MixRound>
    template<typename InputIter>
    void SeedSequence<Count, IntRep, MixRound>::seed(InputIter begin, InputIter end) noexcept {
        assert(std::distance(begin, end) <= Count && "[Error]: Entropy pool size out of range");
        std::copy(begin, end, inits.begin());
        mix_entropy(begin, end);

        // Do additional mixing if and only if for size is very small
        for (size_t i = 1; i < MixRound; ++i)
            mix_entropy(mixer.begin(), mixer.end());
    }

    template<size_t Count, typename IntRep, size_t MixRound>
    template<std::integral T>
    T SeedSequence<Count, IntRep, MixRound>::generate() noexcept {
        if constexpr (std::same_as<T, IntRep>) {
            auto& dataval = mixer[counter];
            dataval ^= hasherB;
            hasherB *= MultB;
            dataval *= hasherB;
            dataval ^= dataval >> XShift;
            counter = (counter + 1) % Count;
            return dataval;
        }
        else {
            static_assert(std::same_as<T, uint64_t>);
            union {
                std::array<uint32_t, 2> i32;
                uint64_t i64;
            } highlow;
            highlow.i32[0] = generate();
            highlow.i32[1] = generate();
            return highlow.i64;
        }
    }

    template<size_t Count, typename IntRep, size_t MixRound>
    template<typename RandomAccessIterator>
    void SeedSequence<Count, IntRep, MixRound>::generate(RandomAccessIterator begin, RandomAccessIterator end) noexcept {
        for (auto it = begin; it != end; ++it)
            *it = generate();
    }

    template<size_t Count, typename IntRep, size_t MixRound>
    template<typename OutputIterator>
    void SeedSequence<Count, IntRep, MixRound>::param(OutputIterator dest) const noexcept {
        std::ranges::copy(inits, dest);
    }

    template<size_t Count, typename IntRep, size_t r>
    void SeedSequence<Count, IntRep, r>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mixer.swap(obj.mixer);
        std::swap(hasherA, obj.hasherA);
        std::swap(hasherB, obj.hasherB);
        std::swap(counter, obj.counter);
    }

    template<size_t Count, typename IntRep, size_t r>
    template<typename InputIter>
    void SeedSequence<Count, IntRep, r>::mix_entropy(InputIter begin, InputIter end) noexcept {
        auto hash = [&](IntRep value) {
            value ^= hasherA;
            hasherA *= MultA;
            value *= hasherA;
            value ^= value >> XShift;
            return value;
        };
        auto mix = [](IntRep x, IntRep y) {
            IntRep result = MixMultL * x - MixMultR * y;
            result ^= result >> XShift;
            return result;
        };

        InputIter current = begin;
        for (auto& elem : mixer) {
            if (current != end)
                elem = hash(*current++);
            else
                elem = hash(0U);
        }

        for (auto& src : mixer)
            for (auto& dest : mixer)
                if (&src != &dest)
                    dest = mix(dest, hash(src));

        for (; current != end; ++current)
            for (auto& dest : mixer)
                dest = mix(dest, hash(*current));
    }
}
/**
 * Random-Number Utilities (randutil)
 *     Addresses common issues with C++11 random number generation.
 *     Makes good seeding easier, and makes using RNGs easy while retaining
 *     all the power.
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2015-2022 Melissa E. O'Neill
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
