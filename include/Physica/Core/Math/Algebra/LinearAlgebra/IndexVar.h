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
    concept IndexVar = std::integral<T> || std::same_as<T, Var> || std::same_as<T, decltype(var())>;

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
            : isFree{std::integral<Ts>...}
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
