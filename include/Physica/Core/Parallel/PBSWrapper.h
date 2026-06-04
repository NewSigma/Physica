/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    /**
     * Target on OpenPBS(https://github.com/openpbs/openpbs)
     */
    class PHYSICA_API PBSWrapper {
        constexpr static size_t hostLength = 64; //64 is enough to hold ipv6

        unsigned int jobCore;
        Array<std::string> hostList;
    public:
        PBSWrapper(const PBSWrapper&) = delete;
        PBSWrapper(PBSWrapper&&) noexcept = delete;
        ~PBSWrapper() = default;
        /* Operators */
        PBSWrapper& operator=(const PBSWrapper&) = delete;
        PBSWrapper& operator=(PBSWrapper&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] unsigned int getJobCore() const noexcept { return jobCore; }
        [[nodiscard]] const Array<std::string>& getHostList() const noexcept { return hostList; }
        /* Static members */
        static const PBSWrapper& getInstance() noexcept; // No [[nodiscard]] for initialization in multi-thread mode
    private:
        PBSWrapper();

        void readJobCore();
        void readHostList();
    };
}
