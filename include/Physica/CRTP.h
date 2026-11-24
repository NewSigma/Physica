/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Macro.h"

namespace Physica {
    /**
     * This class helps implementing CRTP(Curiously Recurring Template Pattern).
     * Typically, classes derived from \class CRTP are named with a 'Mixin' suffix", for example \class ScalarMixin, \class LayerMixin.
     */
    template<class T>
    class CRTP {
        using U = Traits<T>::Derived;
    public:
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr U& getDerived() noexcept { return *static_cast<U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr const U& getDerived() const noexcept { return *static_cast<const U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr U& getConstCastDerived() const noexcept { return *static_cast<U*>(const_cast<CRTP*>(this)); }
    protected:
        constexpr CRTP() = default;
        constexpr CRTP(const CRTP&) = default;
        constexpr CRTP(CRTP&&) noexcept = default;
        ~CRTP() = default;
        /* Operators */
        constexpr CRTP& operator=(const CRTP&) = default;
        constexpr CRTP& operator=(CRTP&&) noexcept = default;
    };
}
