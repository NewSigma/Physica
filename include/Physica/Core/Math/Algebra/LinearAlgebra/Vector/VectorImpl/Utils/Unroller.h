/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    /**
     * Abstraction for a set of logical SIMD registers
     */
    template<Packet Pack, int NumUnroll>
    class [[nodiscard]] Unroller {
        using This = Unroller<Pack, NumUnroll>;
        static_assert(NumUnroll > 0 && (NumUnroll & (NumUnroll - 1)) == 0, "[Error]: Must be power of 2");

        Array<Pack, NumUnroll> packs;
    public:
        Unroller();
        explicit Unroller(std::invocable<size_t> auto generator);
        Unroller(const This&) = default;
        Unroller(This&&) = default;
        ~Unroller() = default;
        /* Operators */
        This& operator=(const This& v) = default;
        This& operator=(This&& v) noexcept = default;
        /* Operations */
        [[nodiscard]] size_t loop_recursive(std::invocable<Pack&, size_t> auto body, size_t length) noexcept;
        [[nodiscard]] Pack sum() noexcept;
    private:
        [[nodiscard]] static size_t loop_recursive(Pack* packs, std::invocable<Pack&, size_t> auto body, size_t length, size_t consumed) noexcept;

        template<Packet, int> friend class Unroller;
    };

    template<Packet Pack, int NumUnroll>
    Unroller<Pack, NumUnroll>::Unroller() : This([](size_t) noexcept { return Pack::zeros(); }) {}

    template<Packet Pack, int NumUnroll>
    Unroller<Pack, NumUnroll>::Unroller(std::invocable<size_t> auto generator)
            : packs(Array<Pack, NumUnroll>::generate(std::move(generator))) {}

    template<Packet Pack, int NumUnroll>
    size_t Unroller<Pack, NumUnroll>::loop_recursive(std::invocable<Pack&, size_t> auto body, size_t length) noexcept {
        return loop_recursive(packs.data(), std::move(body), length, 0);
    }

    template<Packet Pack, int NumUnroll>
    Pack Unroller<Pack, NumUnroll>::sum() noexcept {
        int i = NumUnroll;
        while (i > 1) {
            i /= 2;
            for (int j = 0; j < i; ++j)
                packs[j] += packs[i + j];
        }
        return packs[0];
    }

    template<Packet Pack, int NumUnroll>
    size_t Unroller<Pack, NumUnroll>::loop_recursive(Pack* packs, std::invocable<Pack&, size_t> auto body, const size_t length, size_t consumed) noexcept {
        constexpr size_t StepSize = Pack::size() * NumUnroll;
        for (; consumed < length / StepSize * StepSize; consumed += StepSize) {
            #pragma unroll
            for (int packid = 0; packid < NumUnroll; ++packid)
                packs[packid] = body(packs[packid], consumed + Pack::size() * packid);
        }

        if constexpr (NumUnroll == 1)
            return consumed;
        else
            return Unroller<Pack, NumUnroll / 2>::loop_recursive(packs, std::move(body), length, consumed);
    }
}
