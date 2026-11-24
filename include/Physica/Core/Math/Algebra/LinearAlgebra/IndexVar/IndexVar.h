/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <array>
#include <concepts>
#include "Physica/Macro.h"

namespace Physica {
    /**
     * Dummy class for abstract index notation
     */
    struct [[nodiscard]] Var {
        Var() = default;
        Var(const Var&) = default;
        Var(Var&&) noexcept = default;
        ~Var() = default;
        /* Operators */
        Var& operator=(const Var&) = default;
        Var& operator=(Var&&) noexcept = default;
    };

    [[nodiscard, gnu::nodebug]] __host__ __device__ constexpr auto var() noexcept {
        struct Anonymous : public Var {
            auto operator&() const = delete;
        };
        return Anonymous{};
    }

    template<class T>
    concept IndexVar = std::same_as<T, size_t> || std::same_as<T, Var> || std::same_as<T, decltype(var())>;

    template<IndexVar... Ts>
    struct IndexVarInfo {
        using This = IndexVarInfo;
        constexpr static int NumIndex = sizeof...(Ts);

        std::array<bool, NumIndex> isFree;
        std::array<bool, NumIndex> isNamed;
        std::array<bool, NumIndex> isAnonymous;

        constexpr IndexVarInfo() noexcept;
        constexpr IndexVarInfo(const This&) = default;
        constexpr IndexVarInfo(This&&) noexcept = default;
        constexpr ~IndexVarInfo() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Getters */
        [[nodiscard]] constexpr int getNumFree() const noexcept;
        [[nodiscard]] constexpr int getNumNamed() const noexcept;
        [[nodiscard]] constexpr int getNumAnonymous() const noexcept;
        /* Static members */
        [[nodiscard]] constexpr static int size() noexcept { return NumIndex; }
    };

    template<IndexVar... Ts>
    constexpr IndexVarInfo<Ts...>::IndexVarInfo() noexcept
            : isFree{std::same_as<Ts, size_t>...}
            , isNamed{std::same_as<Ts, Var>...}
            , isAnonymous{std::same_as<Ts, decltype(var())>...} {
        assert((getNumNamed() == 0 || getNumAnonymous() == 0) && "[Error]: It is undefined to mix named and anonymous vars");
    }

    template<IndexVar... Ts>
    constexpr int IndexVarInfo<Ts...>::getNumFree() const noexcept {
        int count = 0;
        for (bool elem : isFree)
            count += elem;
        return count;
    }

    template<IndexVar... Ts>
    constexpr int IndexVarInfo<Ts...>::getNumNamed() const noexcept {
        int count = 0;
        for (bool elem : isNamed)
            count += elem;
        return count;
    }

    template<IndexVar... Ts>
    constexpr int IndexVarInfo<Ts...>::getNumAnonymous() const noexcept {
        int count = 0;
        for (bool elem : isAnonymous)
            count += elem;
        return count;
    }
}
