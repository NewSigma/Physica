/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <string>
#include "Physica/Core/Utils/Container/RingBuffer.h"

namespace Physica {
    /**
     * LogThread will scan and delete buffers that should be deleted.
     */
    class LogBuffer : public RingBuffer {
        bool shouldDelete = false;
    public:
        LogBuffer(size_t size) : RingBuffer(size) {}
        /* Operations */
        std::string makeMsgString();
        std::string formatToString(const char* __restrict format);
        size_t getMsgSize(const char* __restrict format) const;
        void schedualDelete() noexcept { shouldDelete = true; }
        /* Getters */
        [[nodiscard]] bool getShouldDelete() const noexcept { return shouldDelete; }
    };
}
