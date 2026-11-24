/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <cuda_runtime_api.h>

int main() {
    int x = 0, y = 0;
    cudaDeviceGetAttribute(&x, cudaDevAttrConcurrentManagedAccess, 0);
    cudaDeviceGetAttribute(&y, cudaDevAttrPageableMemoryAccess, 0);
    return x & y;
}
