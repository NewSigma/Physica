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
#include <Physica/Python/CXXType.h>
#include <Physica/Python/PhysicaPython.h>

namespace Physica::Python {
    CXXType::CXXType(clang::CXXRecordDecl* pDecl_) : pDecl(pDecl_) {
        assert(pDecl != nullptr && "[Error]: Invalid param");
        const auto& pp = PhysicaPython::getInstance();
        const auto& ctx = pp.getClang().getASTContext();
        const auto type = ctx.getRecordType(pDecl);
        ffiType.size = ctx.getTypeSize(type);
        ffiType.alignment = ctx.getTypeAlign(type);
    }

    CXXType::CXXType(ffi_type ffiType_) : pDecl(nullptr), ffiType(ffiType_) {}

    py::object CXXType::toPython(Core::PlainPtr&& data) const {
        void* pData = data.get();
        switch (ffiType.type) {
        case FFI_TYPE_VOID:
            return py::none();
        case FFI_TYPE_FLOAT:
            return py::float_(*reinterpret_cast<float*>(pData));
        case FFI_TYPE_DOUBLE:
            return py::float_(*reinterpret_cast<double*>(pData));
        default:
            throw std::runtime_error("[Error]: Unknown type");
        }
    }

    Core::PlainPtr CXXType::allocate() const noexcept {
        return Core::PlainPtr(std::aligned_alloc(getAlign(), getSize()));
    }

    void CXXType::swap(CXXType& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pDecl, obj.pDecl);
        std::swap(ffiType, obj.ffiType);
    }
}
