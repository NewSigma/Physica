/*
 * Copyright 2024 Weibo He.
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

#include <string>
#include <utility>

namespace Physica {
    class CXXPtr {
        using This = CXXPtr;
        void* p;
    public:
        CXXPtr(void* p_) : p(p_) {}
        CXXPtr(const This&) = default;
        CXXPtr(This&&) noexcept = default;
        ~CXXPtr() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(void* p1) const noexcept { return p == p1; }
        [[nodiscard]] bool operator!=(void* p1) const noexcept { return p != p1; }
        friend std::ostream& operator<<(std::ostream& os, const This& obj) { return os << obj.p; }
        /* Operations */
        [[nodiscard]] std::string toString() const;
        void swap(CXXPtr& obj) noexcept { std::swap(p, obj.p); }
        /* Getters */
        [[nodiscard]] void* get() const noexcept { return p; }
        template<class T>
        [[nodiscard]] T* cast() noexcept { return reinterpret_cast<T*>(p); }
        template<class T>
        [[nodiscard]] const T* cast() const noexcept { return const_cast<This&>(*this).cast<T>(); }
    };
}
