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
#include <cassert>
#include <Physica/Python/FFI/FuncInfo.h>
#include <Physica/Python/Exception/FFIException.h>

namespace Physica::Python {
    FuncInfo::FuncInfo(unsigned int nargs, const ffi_type* rtype, const ffi_type** atypes) {
        check(ffi_prep_cif(&raw_cif, FFI_DEFAULT_ABI, nargs, const_cast<ffi_type*>(rtype), const_cast<ffi_type**>(atypes)));
    }

    void FuncInfo::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(raw_cif, obj.raw_cif);
    }
}
