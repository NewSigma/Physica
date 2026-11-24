/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <cuda_runtime_api.h>

int main() {
    int result = 0;
    cudaDeviceGetAttribute(&result, cudaDevAttrWarpSize, 0);
    return result;
}
