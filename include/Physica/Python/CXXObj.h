/*
 * Copyright 2024 WeiBo He.
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

#include <utility>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace clang {
    class CXXRecordDecl;
}

namespace Physica::Python {
    class CXXObj {
        using This = CXXObj;
    private:
        const clang::CXXRecordDecl* decl;
        void* pObj;
    public:
        CXXObj() = default;
        CXXObj(const CXXObj&) = delete;
        CXXObj(CXXObj&&) noexcept;
        ~CXXObj();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        py::object call(py::handle type, const char* name);

        template<class T>
        [[nodiscard]] T& getDerived() noexcept { return *reinterpret_cast<T*>(pObj); }
        template<class T>
        [[nodiscard]] const T& getDerived() const noexcept { return const_cast<This&>(*this).getDerived<T>(); }
        void swap(CXXObj& __restrict obj) noexcept;
    private:
        clang::FunctionDecl* makeDestructorAST();
    };

    inline void CXXObj::swap(CXXObj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(decl, obj.decl);
        std::swap(pObj, obj.pObj);
    }
}
