/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Exception/MKL/VSL.h"

using namespace Physica;

namespace {
    class Impl final : public std::error_category {
    public:
        Impl() = default;
        Impl(const Impl&) = delete;
        Impl(Impl&&) noexcept = delete;
        ~Impl() = default;
        /* Operators */
        Impl& operator=(const Impl&) = delete;
        Impl& operator=(Impl&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] const char* name() const noexcept final { return "Intel MKL VSL"; }
        [[nodiscard]] std::string message(int) const final;
    };

    std::string Impl::message(int code) const {
        switch (code) {
        case VSL_STATUS_OK:
            return "No error";
        case VSL_ERROR_FEATURE_NOT_IMPLEMENTED:
            return "Feature not implemented";
        case VSL_ERROR_UNKNOWN:
            return "Unknown internal error";
        case VSL_ERROR_BADARGS:
            return "Invalid arguments";
        case VSL_ERROR_MEM_FAILURE:
            return "Memory allocation failure";
        case VSL_ERROR_NULL_PTR:
            return "Null pointer";
        case VSL_ERROR_CPU_NOT_SUPPORTED:
            return "CPU not supported";
        case VSL_RNG_ERROR_INVALID_BRNG_INDEX:
            return "Invalid BRNG index";
        case VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED:
            return "Leapfrog method not supported by the BRNG";
        case VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED:
            return "Skipahead method not supported by the BRNG";
        case VSL_RNG_ERROR_SKIPAHEADEX_UNSUPPORTED:
            return "SkipaheadEx method not supported by the BRNG";
        case VSL_RNG_ERROR_BRNGS_INCOMPATIBLE:
            return "BRNGs are incompatible";
        case VSL_RNG_ERROR_BAD_STREAM:
            return "Invalid stream pointer";
        case VSL_RNG_ERROR_BRNG_TABLE_FULL:
            return "BRNG table is full";
        case VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE:
            return "Invalid stream state size";
        case VSL_RNG_ERROR_BAD_WORD_SIZE:
            return "Invalid word size";
        case VSL_RNG_ERROR_BAD_NSEEDS:
            return "Invalid number of seeds";
        case VSL_RNG_ERROR_BAD_NBITS:
            return "Invalid number of bits";
        case VSL_RNG_ERROR_QRNG_PERIOD_ELAPSED:
            return "QRNG period elapsed";
        case VSL_RNG_ERROR_LEAPFROG_NSTREAMS_TOO_BIG:
            return "Number of streams in leapfrog is too big";
        case VSL_RNG_ERROR_BRNG_NOT_SUPPORTED:
            return "BRNG not supported";
        case VSL_RNG_ERROR_BAD_UPDATE:
            return "Bad update";
        case VSL_RNG_ERROR_NO_NUMBERS:
            return "No numbers generated";
        case VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM:
            return "Invalid abstract stream";
        case VSL_RNG_ERROR_NONDETERM_NOT_SUPPORTED:
            return "Non-deterministic generator not supported";
        case VSL_RNG_ERROR_NONDETERM_NRETRIES_EXCEEDED:
            return "Number of retries for non-deterministic generator exceeded";
        case VSL_RNG_ERROR_ARS5_NOT_SUPPORTED:
            return "ARS5 generator not supported";
        case VSL_DISTR_MULTINOMIAL_BAD_PROBABILITY_ARRAY:
            return "Category probabilities for multinomial distribution do not sum to 1.0";
        case VSL_RNG_ERROR_FILE_CLOSE:
            return "File close error";
        case VSL_RNG_ERROR_FILE_OPEN:
            return "File open error";
        case VSL_RNG_ERROR_FILE_WRITE:
            return "File write error";
        case VSL_RNG_ERROR_FILE_READ:
            return "File read error";
        case VSL_RNG_ERROR_BAD_FILE_FORMAT:
            return "Bad file format";
        case VSL_RNG_ERROR_UNSUPPORTED_FILE_VER:
            return "Unsupported file version";
        case VSL_RNG_ERROR_BAD_MEM_FORMAT:
            return "Bad memory format";
        default:
            return "Bad error code";
        }
    }
}

VSLException::VSLException(int code) noexcept : std::system_error(code, Impl()) {}
