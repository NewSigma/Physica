/*
 * Copyright 2021-2026 Weibo He.
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
module;

#include <string>
#include "Physica/Core/Utils/Container/RingBuffer.h"

export module Physica.Logger.LogBuffer;

export namespace Physica {
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
