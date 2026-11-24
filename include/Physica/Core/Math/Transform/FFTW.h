/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include "Physica/Macro.h"
#include "Physica/Core/Utils/Handle.h"

namespace Physica {
    using FFTPlan = Handle<HandleType::FFTW_Plan>;

    enum class PlanFlag : int8_t {
        Measure,
        Estimate,
        Patient,
        Exhaustive
    };

    namespace Internal {
        PHYSICA_API void* fftw_malloc(size_t size) noexcept;
        PHYSICA_API void fftw_free(void* p) noexcept;

        template<bool isSinglePrec>
        [[nodiscard]] PHYSICA_API FFTPlan fftw_plan_dft(int rank, const int* dims, void* buffer, bool forward, PlanFlag flag) noexcept;
        template<bool isSinglePrec>
        [[nodiscard]] PHYSICA_API FFTPlan fftw_plan_dft_r2c(int rank, const int* dims, void* buffer, PlanFlag flag) noexcept;
        template<bool isSinglePrec>
        [[nodiscard]] PHYSICA_API FFTPlan fftw_plan_dft_c2r(int rank, const int* dims, void* buffer, PlanFlag flag) noexcept;

        template<bool isSinglePrec>
        PHYSICA_API void fftw_execute_dft(FFTPlan plan, void* buffer) noexcept;
        template<bool isSinglePrec>
        PHYSICA_API void fftw_execute_dft_r2c(FFTPlan plan, void* buffer) noexcept;
        template<bool isSinglePrec>
        PHYSICA_API void fftw_execute_dft_c2r(FFTPlan plan, void* buffer) noexcept;
        template<bool isSinglePrec>
        PHYSICA_API void fftw_destroy_plan(FFTPlan plan) noexcept;

        extern template FFTPlan fftw_plan_dft<false>(int, const int*, void*, bool, PlanFlag) noexcept;
        extern template FFTPlan fftw_plan_dft<true>(int, const int*, void*, bool, PlanFlag) noexcept;
        extern template FFTPlan fftw_plan_dft_r2c<false>(int, const int*, void*, PlanFlag) noexcept;
        extern template FFTPlan fftw_plan_dft_r2c<true>(int, const int*, void*, PlanFlag) noexcept;
        extern template FFTPlan fftw_plan_dft_c2r<false>(int, const int*, void*, PlanFlag) noexcept;
        extern template FFTPlan fftw_plan_dft_c2r<true>(int, const int*, void*, PlanFlag) noexcept;

        extern template void fftw_execute_dft<false>(FFTPlan, void*) noexcept;
        extern template void fftw_execute_dft<true>(FFTPlan, void*) noexcept;
        extern template void fftw_execute_dft_r2c<false>(FFTPlan, void*) noexcept;
        extern template void fftw_execute_dft_r2c<true>(FFTPlan, void*) noexcept;
        extern template void fftw_execute_dft_c2r<false>(FFTPlan, void*) noexcept;
        extern template void fftw_execute_dft_c2r<true>(FFTPlan, void*) noexcept;
        extern template void fftw_destroy_plan<false>(FFTPlan) noexcept;
        extern template void fftw_destroy_plan<true>(FFTPlan) noexcept;
    }
}
