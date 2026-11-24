/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cassert>
#include <cstddef>
#include <memory>
#include <ranges>
#include "H5Type.h"

namespace Physica {
    class PHYSICA_API H5Attribute : public H5ID {
        using This = H5Attribute;
    public:
        H5Attribute() = default;
        explicit H5Attribute(H5ID id) noexcept;
        H5Attribute(const This&) = default;
        H5Attribute(This&&) noexcept = default;
        ~H5Attribute() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void read(auto& value) const;
        void write(const auto& value) const;

        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept;
        /* Static members */
        [[nodiscard]] constexpr static IdentifierType itype() noexcept { return IdentifierType::Attribute; }
    private:
        void read(const H5Type& dtype, void* buf) const;
        void write(const H5Type& dtype, const void* buf) const;
    };

    void H5Attribute::read(auto& value) const {
        using T = std::remove_cvref_t<decltype(value)>;
        if constexpr (std::ranges::contiguous_range<T>) {
            assert(std::ranges::size(value) == getSize() && "[Error]: Inconsistent range size");
            read(H5Type::get<T>(), std::ranges::data(value));
        }
        else
            read(H5Type::get<T>(), std::addressof(value));
    }

    void H5Attribute::write(const auto& value) const {
        using T = std::remove_cvref_t<decltype(value)>;
        if constexpr (std::ranges::contiguous_range<T>) {
            assert(std::ranges::size(value) == getSize() && "[Error]: Inconsistent range size");
            write(H5Type::get<T>(), std::ranges::data(value));
        }
        else
            write(H5Type::get<T>(), std::addressof(value));
    }
}
