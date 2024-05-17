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

#include <type_traits>
#include <typeinfo>
#include <utility>

namespace Physica::Python {
    class CXXObj {
        using This = CXXObj;

        const char* typeName;
        void* pObj;
    public:
        template<class T>
        CXXObj(T&& obj) noexcept;
        CXXObj(const CXXObj&) = delete;
        CXXObj(CXXObj&&) noexcept;
        ~CXXObj();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class T>
        [[nodiscard]] T& getDerived() noexcept { return *reinterpret_cast<T*>(pObj); }
        template<class T>
        [[nodiscard]] const T& getDerived() const noexcept { return const_cast<This&>(*this).getDerived<T>(); }
        /* Getter */
        [[nodiscard]] const char* getTypeName() const noexcept { return typeName; }
    };

    template<class T>
    CXXObj::CXXObj(T&& obj) noexcept {
        using T1 = typename std::remove_cv<T>::type;
        using T2 = typename std::remove_reference<T1>::type;
        typeName = typeid(T2).name();
        pObj = reinterpret_cast<void*>(new T(std::move(obj)));
    }
    
    template<class T>
    inline CXXObj toPython(T&& obj) noexcept {
        return CXXObj(std::move(obj));
    }
}
