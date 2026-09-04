/*
 * Copyright 2023-2026 Weibo He.
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

#include "../CompactMatrix.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    class Flatten<M> : public CompactVector<Flatten<M>> {
        using This = Flatten<M>;

        decay_rvalue_t<M> mat;
    public:
        using Base = CompactVector<This>;
    public:
        Flatten(M&& mat_) : mat(std::forward<M>(mat_)) {}
        Flatten(const This&) = default;
        Flatten(This&&) noexcept = default;
        ~Flatten() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) grads(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getSize(); }
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    decltype(auto) Flatten<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().flatten();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    decltype(auto) Flatten<M>::grads(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).grads().flatten();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    auto Flatten<M>::data(this auto&& self) noexcept {
        return self.mat.data();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ consteval size_t Flatten<M>::getSizeAtCompile() noexcept {
        return std::remove_reference_t<M>::getSizeAtCompile();
    }
}
