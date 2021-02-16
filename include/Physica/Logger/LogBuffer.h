/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Utils/RingBuffer.h"

namespace Physica::Logger {
    /**
     * LogThread will scan and delete buffers that should be deleted.
     */
    class LogBuffer : public Utils::RingBuffer {
        bool shouldDelete;
    public:
        LogBuffer(size_t size) : RingBuffer(size), shouldDelete(false) {}

        void schedualDelete() noexcept { shouldDelete = true; }
        [[nodiscard]] bool getShouldDelete() const noexcept { return shouldDelete; }
    };
}
