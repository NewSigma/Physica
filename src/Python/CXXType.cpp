/*
 * Copyright 2024-2025 Weibo He.
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
#include <clang/AST/RecordLayout.h>
#include "Physica/Python/CXXType.h"
#include "Physica/Python/PhysicaPython.h"

using namespace Physica;

CXXType::CXXType(clang::CXXRecordDecl* pDecl_) : pDecl(pDecl_) {
    assert(pDecl != nullptr && "[Error]: Invalid param");
    auto& pp = PhysicaPython::getInstance();
    const auto& ctx = pp.getClang().getASTContext();
    const auto& layout = ctx.getASTRecordLayout(pDecl);
    ffiType.size = layout.getSize().getQuantity();
    ffiType.alignment = layout.getAlignment().getQuantity();
}

CXXType::CXXType(ffi_type ffiType_)
        : pDecl(nullptr), ffiType(ffiType_) {}

nanobind::object CXXType::toPython(void* data) const {
    switch (ffiType.type) {
    case FFI_TYPE_VOID:
        return nanobind::none();
    case FFI_TYPE_FLOAT:
        return nanobind::float_(*reinterpret_cast<float*>(data));
    case FFI_TYPE_DOUBLE:
        return nanobind::float_(*reinterpret_cast<double*>(data));
    default:
        throw std::runtime_error("[Error]: Unknown type");
    }
}

void CXXType::swap(CXXType& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(pDecl, obj.pDecl);
    std::swap(ffiType, obj.ffiType);
}
