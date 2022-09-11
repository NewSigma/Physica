/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    /**
     * Target on OpenPBS(https://github.com/openpbs/openpbs)
     */
    class PBSWarpper {
        constexpr static size_t hostLength = 64; //64 is enough to hold ipv6

        unsigned int jobCore;
        Utils::Array<std::string> hostList;
    public:
        PBSWarpper(const PBSWarpper&) = delete;
        PBSWarpper(PBSWarpper&&) noexcept = delete;
        ~PBSWarpper() = default;
        /* Operators */
        PBSWarpper& operator=(const PBSWarpper&) = delete;
        PBSWarpper& operator=(PBSWarpper&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] unsigned int getJobCore() const noexcept { return jobCore; }
        [[nodiscard]] const Utils::Array<std::string>& getHostList() const noexcept { return hostList; }
        /* Static members */
        [[nodiscard]] static const PBSWarpper& getInstance();
    private:
        PBSWarpper();

        void readJobCore();
        void readHostList();
    };
}
