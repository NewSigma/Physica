/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
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
