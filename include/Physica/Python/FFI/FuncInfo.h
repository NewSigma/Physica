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

#include <ffi.h>

namespace Physica::Python {
    class FuncInfo {
        using This = FuncInfo;

        ffi_cif raw_cif;
    public:
        FuncInfo(unsigned int nargs, const ffi_type* rtype, const ffi_type** atypes);
        FuncInfo(const FuncInfo&) = default;
        FuncInfo(FuncInfo&&) = default;
        ~FuncInfo() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const ffi_cif* cif() const noexcept { return &raw_cif; }
    };
}
