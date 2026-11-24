/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/IO/HDF5/H5Type.h"
#include <cassert>
#include "Physica/Core/Utils/Builtin.h"

using namespace Physica;

H5Type::H5Type(H5ID id) noexcept : H5ID(std::move(id)) {
    assert(isa<H5Type>());
}

H5Type::H5Type(hid_t hid) noexcept : H5Type(H5ID(hid)) {}

void H5Type::insert(const char* name, size_t offset, const H5Type& memberType) noexcept {
    assert(isCompound());
    H5Tinsert(getHID(), name, offset, memberType.getHID());
}

bool H5Type::isCompound() const noexcept {
    return H5Tget_class(getHID()) == H5T_COMPOUND;
}

auto H5Type::getPrimitiveType(PrimitiveType type) noexcept -> This {
    switch (type) {
    case PrimitiveType::Int8:
        return This(H5T_NATIVE_INT8);
    case PrimitiveType::Int16:
        return This(H5T_NATIVE_INT16);
    case PrimitiveType::Int32:
        return This(H5T_NATIVE_INT32);
    case PrimitiveType::Int64:
        return This(H5T_NATIVE_INT64);
    case PrimitiveType::UInt8:
        return This(H5T_NATIVE_UINT8);
    case PrimitiveType::UInt16:
        return This(H5T_NATIVE_UINT16);
    case PrimitiveType::UInt32:
        return This(H5T_NATIVE_UINT32);
    case PrimitiveType::UInt64:
        return This(H5T_NATIVE_UINT64);
    case PrimitiveType::Bool:
        return This(H5T_NATIVE_HBOOL);
    case PrimitiveType::Char:
        return This(H5T_NATIVE_CHAR);
    case PrimitiveType::SignedChar:
        return This(H5T_NATIVE_SCHAR);
    case PrimitiveType::UnsignedChar:
        return This(H5T_NATIVE_UCHAR);
    case PrimitiveType::Short:
        return This(H5T_NATIVE_SHORT);
    case PrimitiveType::UnsignedShort:
        return This(H5T_NATIVE_USHORT);
    case PrimitiveType::Int:
        return This(H5T_NATIVE_INT);
    case PrimitiveType::UnsignedInt:
        return This(H5T_NATIVE_UINT);
    case PrimitiveType::Long:
        return This(H5T_NATIVE_LONG);
    case PrimitiveType::UnsignedLong:
        return This(H5T_NATIVE_ULONG);
    case PrimitiveType::LongLong:
        return This(H5T_NATIVE_LLONG);
    case PrimitiveType::UnsignedLongLong:
        return This(H5T_NATIVE_ULLONG);
    case PrimitiveType::Float:
        return This(H5T_NATIVE_FLOAT);
    case PrimitiveType::Double:
        return This(H5T_NATIVE_DOUBLE);
    case PrimitiveType::LongDouble:
        return This(H5T_NATIVE_LDOUBLE);
    default:
        unreachable();
    }
}
