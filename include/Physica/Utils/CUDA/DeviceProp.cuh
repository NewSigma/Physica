/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Macro.h>
#include <Physica/Utils/Container/Array/Array.h>

namespace Physica::Utils {
    class PHYSICA_API DeviceProp {
    private:
        int driverVersion;
        int runtimeVersion;
        Array<cudaDeviceProp> propList;
    public:
        DeviceProp(const DeviceProp&) = delete;
        DeviceProp(DeviceProp&&) noexcept = delete;
        ~DeviceProp() = default;
        /* Operators */
        DeviceProp& operator=(const DeviceProp&) = delete;
        DeviceProp& operator=(DeviceProp&&) noexcept = delete;
        /* Operations */
        std::ostream& printDeviceProp(unsigned int device, std::ostream& os) const;
        /* Getters */
        [[nodiscard]] int getDriverVersion() const { return driverVersion; }
        [[nodiscard]] int getRuntimeVersion() const { return runtimeVersion; }
        [[nodiscard]] const cudaDeviceProp& getProperty(size_t device) const { return propList[device]; }
        /* Static members */
        static const DeviceProp& getInstance(); // No [[nodiscard]] for initialization in multi-thread mode
    private:
        DeviceProp();
    };
}
