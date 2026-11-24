/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cuda_runtime_api.h>

int main() {
    int driverVersion = 0;
    int deviceSupportsMemoryPools = 0;
    cudaDriverGetVersion(&driverVersion);
    if (driverVersion >= 11020)
        cudaDeviceGetAttribute(&deviceSupportsMemoryPools, cudaDevAttrMemoryPoolsSupported, 0);
    return deviceSupportsMemoryPools;
}
