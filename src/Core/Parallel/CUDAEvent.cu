/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <cassert>
#include "Physica/Core/Parallel/CUDAEvent.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"

using namespace Physica;

CUDAEvent::CUDAEvent(int flag) {
    check(cudaEventCreateWithFlags(&event, flag));
}

CUDAEvent::CUDAEvent(CUDAEvent&& obj) noexcept : event(obj.event) {
    obj.event = nullptr;
}

void CUDAEvent::record() {
    cudaEventRecord(event, CUDAContext::getInstance());
}

CUDAEvent::~CUDAEvent() {
    cudaEventDestroy(event);
}

void CUDAEvent::swap(CUDAEvent& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(event, obj.event);
}
