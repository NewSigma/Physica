/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) grads(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return mat.getSize(); }
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
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
    auto Flatten<M>::data_handle(this auto&& self) noexcept {
        return self.mat.data_handle();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ consteval size_t Flatten<M>::getSizeAtCompile() noexcept {
        return std::remove_reference_t<M>::getSizeAtCompile();
    }
}
