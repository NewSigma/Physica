/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <limits>
#include <stack>
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"
#include "Physica/Core/Exception/CUDA/cuDSS.cuh"
#include "Physica/Core/Exception/CUDA/cuSolver.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"

using namespace Physica;

CUDAContext::Page::Page(int device_, CUDAStream stream_) : device(device_), stream(std::move(stream_)) {
    check(cublasCreate(&cublas));
    check(cublasSetStream(cublas, stream));

    check(cusolverDnCreate(&cuSolverDn));
    check(cusolverDnSetStream(cuSolverDn, stream));
    check(cusolverDnCreateParams(&cuSolverDnParams));

    check(cudssCreate(&cudss));
    check(cudssSetStream(cudss, stream));
}

CUDAContext::Page::~Page() {
    cublasDestroy(cublas);
    cusolverDnDestroy(cuSolverDn);
    cusolverDnDestroyParams(cuSolverDnParams);
    cudssDestroy(cudss);
}

CUDAContext::CUDAContext() {
    setDevice(0);
    pages.emplace(0, CUDAStream(nullptr));
}

CUDAContext::~CUDAContext() {
    while (!pages.empty())
        pages.pop();
}

CUDAContext::operator cudaStream_t() const noexcept {
    return stream();
}

CUDAContext::operator cublasContext*() const noexcept {
    return pages.top().cublas;
}

CUDAContext::operator cusolverDnContext*() const noexcept {
    return pages.top().cuSolverDn;
}

CUDAContext::operator cusolverDnParams*() const noexcept {
    return pages.top().cuSolverDnParams;
}

CUDAContext::operator cudssContext*() const noexcept {
    return pages.top().cudss;
}

auto CUDAContext::push() -> PageGuard {
    return push(device());
}

auto CUDAContext::push(int device) -> PageGuard {
    setDevice(0);
    pages.emplace(device, CUDAStream());
    return {};
}

void CUDAContext::wait() const {
    stream().wait();
}

int CUDAContext::device() const noexcept {
    return pages.top().device;
}

const CUDAStream& CUDAContext::stream() const noexcept {
    return pages.top().stream;
}

CUDAContext& CUDAContext::getInstance() noexcept {
    thread_local static CUDAContext instance{};
    return instance;
}

void CUDAContext::setPointerMode(bool isDeviceSide) const noexcept {
    cublasSetPointerMode(*this, isDeviceSide ? CUBLAS_POINTER_MODE_DEVICE : CUBLAS_POINTER_MODE_HOST);
}

void CUDAContext::setDevice(int device) {
    check(cudaSetDevice(device));

    uint64_t val = std::numeric_limits<uint64_t>::max();
    cudaMemPool_t memPool{};
    check(cudaDeviceGetDefaultMemPool(&memPool, device));
    check(cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &val));
}
