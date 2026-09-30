/*
 * Copyright 2026 Weibo He.
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
#include <bit>
#include <fftw3.h>
#include "Physica/Core/Math/Transform/FFTW.h"

using namespace Physica;

namespace {
    [[nodiscard]] unsigned int toFlag(PlanFlag flag) noexcept {
        switch (flag) {
        case PlanFlag::Measure:
            return FFTW_MEASURE;
        case PlanFlag::Estimate:
            return FFTW_ESTIMATE;
        case PlanFlag::Patient:
            return FFTW_PATIENT;
        case PlanFlag::Exhaustive:
            return FFTW_EXHAUSTIVE;
        default:
            return FFTW_ESTIMATE;
        }
    }
}

void* Physica::Internal::fftw_malloc(size_t size) noexcept {
    return ::fftw_malloc(size);
}

void Physica::Internal::fftw_free(void* p) noexcept {
    ::fftw_free(p);
}

template<bool isSinglePrec>
FFTPlan Physica::Internal::fftw_plan_dft(int rank, const int* dims, void* buffer, bool forward, PlanFlag flag) noexcept {
    const int sign = forward ? FFTW_FORWARD : FFTW_BACKWARD;
    if constexpr (isSinglePrec)
        return FFTPlan(::fftwf_plan_dft(rank, std::bit_cast<int*>(dims), std::bit_cast<fftwf_complex*>(buffer), std::bit_cast<fftwf_complex*>(buffer), sign, toFlag(flag)));
    else
        return FFTPlan(::fftw_plan_dft(rank, std::bit_cast<int*>(dims), std::bit_cast<fftw_complex*>(buffer), std::bit_cast<fftw_complex*>(buffer), sign, toFlag(flag)));
}

template<bool isSinglePrec>
FFTPlan Physica::Internal::fftw_plan_dft_r2c(int rank, const int* dims, void* buffer, PlanFlag flag) noexcept {
    if constexpr (isSinglePrec)
        return FFTPlan(::fftwf_plan_dft_r2c(rank, std::bit_cast<int*>(dims), std::bit_cast<float*>(buffer), std::bit_cast<fftwf_complex*>(buffer), toFlag(flag)));
    else
        return FFTPlan(::fftw_plan_dft_r2c(rank, std::bit_cast<int*>(dims), std::bit_cast<double*>(buffer), std::bit_cast<fftw_complex*>(buffer), toFlag(flag)));
}

template<bool isSinglePrec>
FFTPlan Physica::Internal::fftw_plan_dft_c2r(int rank, const int* dims, void* buffer, PlanFlag flag) noexcept {
    if constexpr (isSinglePrec)
        return FFTPlan(::fftwf_plan_dft_c2r(rank, std::bit_cast<int*>(dims), std::bit_cast<fftwf_complex*>(buffer), std::bit_cast<float*>(buffer), toFlag(flag)));
    else
        return FFTPlan(::fftw_plan_dft_c2r(rank, std::bit_cast<int*>(dims), std::bit_cast<fftw_complex*>(buffer), std::bit_cast<double*>(buffer), toFlag(flag)));
}

template<bool isSinglePrec>
void Physica::Internal::fftw_execute_dft(FFTPlan plan, void* buffer) noexcept {
    if constexpr (isSinglePrec)
        ::fftwf_execute_dft(std::bit_cast<fftwf_plan>(plan), std::bit_cast<fftwf_complex*>(buffer), std::bit_cast<fftwf_complex*>(buffer));
    else
        ::fftw_execute_dft(std::bit_cast<fftw_plan>(plan), std::bit_cast<fftw_complex*>(buffer), std::bit_cast<fftw_complex*>(buffer));
}

template<bool isSinglePrec>
void Physica::Internal::fftw_execute_dft_r2c(FFTPlan plan, void* buffer) noexcept {
    if constexpr (isSinglePrec)
        ::fftwf_execute_dft_r2c(std::bit_cast<fftwf_plan>(plan), std::bit_cast<float*>(buffer), std::bit_cast<fftwf_complex*>(buffer));
    else
        ::fftw_execute_dft_r2c(std::bit_cast<fftw_plan>(plan), std::bit_cast<double*>(buffer), std::bit_cast<fftw_complex*>(buffer));
}

template<bool isSinglePrec>
void Physica::Internal::fftw_execute_dft_c2r(FFTPlan plan, void* buffer) noexcept {
    if constexpr (isSinglePrec)
        ::fftwf_execute_dft_c2r(std::bit_cast<fftwf_plan>(plan), std::bit_cast<fftwf_complex*>(buffer), std::bit_cast<float*>(buffer));
    else
        ::fftw_execute_dft_c2r(std::bit_cast<fftw_plan>(plan), std::bit_cast<fftw_complex*>(buffer), std::bit_cast<double*>(buffer));
}

template<bool isSinglePrec>
void Physica::Internal::fftw_destroy_plan(FFTPlan plan) noexcept {
    if constexpr (isSinglePrec)
        ::fftwf_destroy_plan(std::bit_cast<fftwf_plan>(plan));
    else
        ::fftw_destroy_plan(std::bit_cast<fftw_plan>(plan));
}

template FFTPlan Physica::Internal::fftw_plan_dft<false>(int, const int*, void*, bool, PlanFlag) noexcept;
template FFTPlan Physica::Internal::fftw_plan_dft<true>(int, const int*, void*, bool, PlanFlag) noexcept;
template FFTPlan Physica::Internal::fftw_plan_dft_r2c<false>(int, const int*, void*, PlanFlag) noexcept;
template FFTPlan Physica::Internal::fftw_plan_dft_r2c<true>(int, const int*, void*, PlanFlag) noexcept;
template FFTPlan Physica::Internal::fftw_plan_dft_c2r<false>(int, const int*, void*, PlanFlag) noexcept;
template FFTPlan Physica::Internal::fftw_plan_dft_c2r<true>(int, const int*, void*, PlanFlag) noexcept;

template void Physica::Internal::fftw_execute_dft<false>(FFTPlan, void*) noexcept;
template void Physica::Internal::fftw_execute_dft<true>(FFTPlan, void*) noexcept;
template void Physica::Internal::fftw_execute_dft_r2c<false>(FFTPlan, void*) noexcept;
template void Physica::Internal::fftw_execute_dft_r2c<true>(FFTPlan, void*) noexcept;
template void Physica::Internal::fftw_execute_dft_c2r<false>(FFTPlan, void*) noexcept;
template void Physica::Internal::fftw_execute_dft_c2r<true>(FFTPlan, void*) noexcept;
template void Physica::Internal::fftw_destroy_plan<false>(FFTPlan) noexcept;
template void Physica::Internal::fftw_destroy_plan<true>(FFTPlan) noexcept;
