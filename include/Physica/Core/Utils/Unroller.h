/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <span>
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

        Array<Pack, NumUnroll> packs{};
    public:
        Unroller() = default;
        explicit Unroller(std::invocable<size_t> auto generator);
        Unroller(const This&) = default;
        Unroller(This&&) noexcept = default;
        ~Unroller() = default;
        /* Operators */
        This& operator=(const This& v) = default;
        This& operator=(This&& v) noexcept = default;
        /* Operations */
        [[nodiscard]] size_t loop_recursive(std::invocable<Pack, size_t> auto body, size_t length) noexcept;
        [[nodiscard]] Pack sum() noexcept;
        /* Static members */
        static Unroller zeros();
    private:
        [[nodiscard]] static size_t loop_recursive(std::span<Pack, NumUnroll> packs, std::invocable<Pack, size_t> auto body, size_t length, size_t consumed) noexcept;

        template<Packet, int> friend class Unroller;
    };

    template<Packet Pack, int NumUnroll>
    Unroller<Pack, NumUnroll>::Unroller(std::invocable<size_t> auto generator)
            : packs(Array<Pack, NumUnroll>::generate(std::move(generator))) {}

    template<Packet Pack, int NumUnroll>
    size_t Unroller<Pack, NumUnroll>::loop_recursive(std::invocable<Pack, size_t> auto body, size_t length) noexcept {
        return loop_recursive(std::span<Pack, NumUnroll>(packs), std::move(body), length, 0);
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
    Unroller<Pack, NumUnroll> Unroller<Pack, NumUnroll>::zeros() {
        return This([](size_t) noexcept { return Pack::zeros(); });
    }

    template<Packet Pack, int NumUnroll>
    size_t Unroller<Pack, NumUnroll>::loop_recursive(std::span<Pack, NumUnroll> packs, std::invocable<Pack, size_t> auto body, const size_t length, size_t consumed) noexcept {
        constexpr size_t StepSize = Pack::size() * NumUnroll;
        for (; consumed < length / StepSize * StepSize; consumed += StepSize) {
            #pragma unroll
            for (int packid = 0; packid < NumUnroll; ++packid)
                packs[packid] = body(packs[packid], consumed + Pack::size() * packid);
        }

        if constexpr (NumUnroll == 1)
            return consumed;
        else
            return Unroller<Pack, NumUnroll / 2>::loop_recursive(packs.template first<NumUnroll / 2>(), std::move(body), length, consumed);
    }
}
