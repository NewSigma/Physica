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

namespace Physica {
    class CXXPtr;

    class CXXObj {
        using This = CXXObj;
        using ForeignFunc = void (*)();
        using CXXRecordDecl = clang::CXXRecordDecl;
        using ClassTemplateDecl = clang::ClassTemplateDecl;

        struct deleter {
            void operator()(void* ptr) { ::operator delete(ptr); }
        };

        using Ptr = std::unique_ptr<void, deleter>;
    private:
        CXXRecordDecl* pDecl;
        Ptr pObj;
    public:
        CXXObj(CXXPtr p, nanobind::args tparams);
        CXXObj(const This&) = delete;
        CXXObj(This&&) noexcept;
        ~CXXObj();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void construct(const nanobind::args& args);
        [[nodiscard]] nanobind::object call(const char* rtnTyName, const char* name, const nanobind::args& args);
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] CXXType getType() const noexcept { return CXXType(pDecl); }
        template<class T>
        [[nodiscard]] auto& getDerived(this auto&& self) noexcept;
    private:
        CXXObj(CXXRecordDecl* pDecl_, void* pObj_) noexcept;

        [[nodiscard]] static CXXRecordDecl* findSpecialization(const ClassTemplateDecl& templateDecl, nanobind::args tparams);
        [[nodiscard]] static bool matchParamT(const nanobind::handle& pyarg, const clang::TemplateArgument& targ);
        [[nodiscard]] static ForeignFunc lookupFunc(clang::GlobalDecl decl);
        [[nodiscard]] static Ptr allocate(CXXType ty);
    };

    inline void CXXObj::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pDecl, obj.pDecl);
        std::swap(pObj, obj.pObj);
    }

    template<class T>
    auto& CXXObj::getDerived(this auto&& self) noexcept {
        if constexpr (std::is_const_v<std::remove_reference_t<decltype(self)>>)
            return *reinterpret_cast<const T*>(self.pObj.get());
        else
            return *reinterpret_cast<T*>(self.pObj.get());
    }
}
