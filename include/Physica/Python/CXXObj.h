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

#include <type_traits>
#include <utility>
#include "PhysicaPython.h"

namespace clang {
    class CXXRecordDecl;
}

namespace Physica::Python {
    class CXXObj {
        using This = CXXObj;
        using CXXRecordDecl = clang::CXXRecordDecl;
    private:
        const CXXRecordDecl* pDecl;
        void* pObj;
    public:
        CXXObj(const CXXObj&) = delete;
        CXXObj(CXXObj&&) noexcept;
        ~CXXObj();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] py::object call(const char* rtnTyName, const char* name, py::args args);
        inline void swap(CXXObj& __restrict obj) noexcept;
        /* Getters */
        template<class T>
        [[nodiscard]] T& getDerived() noexcept { return *reinterpret_cast<T*>(pObj); }
        template<class T>
        [[nodiscard]] const T& getDerived() const noexcept { return const_cast<This&>(*this).getDerived<T>(); }
        /* Static members */
        [[nodiscard]] static This create(const CXXRecordDecl* pDecl);
    private:
        CXXObj(const CXXRecordDecl* pDecl_, void* pObj_) noexcept;

        static void* allocateObj(const CXXRecordDecl* pDecl);
    };

    inline void CXXObj::swap(CXXObj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pDecl, obj.pDecl);
        std::swap(pObj, obj.pObj);
    }
}
