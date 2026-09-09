/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Handle.h"

namespace Physica {
    class PHYSICA_API Request {
        using This = Request;
    public:
        using Handle = Handle<HandleType::MPI_Request>;
    private:
        Handle h;
    public:
        Request();
        Request(Handle h);
        Request(const This&) = delete;
        Request(This&& obj) noexcept;
        ~Request();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void wait();
        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] bool query();
    };
}
