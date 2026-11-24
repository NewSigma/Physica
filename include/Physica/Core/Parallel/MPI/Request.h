/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
