/*
 * Copyright 2024 Weibo He.
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
#include <Physica/Core/Exception/MKL/VSL.h>

#ifdef PHYSICA_MKL

namespace Physica::Core {
    std::string VSLException::Impl::message(int code) const {
        switch(code) {
        case VSL_STATUS_OK:
            return "No error";
        case VSL_ERROR_FEATURE_NOT_IMPLEMENTED:
        case VSL_ERROR_UNKNOWN:
        case VSL_ERROR_BADARGS:
        case VSL_ERROR_MEM_FAILURE:
        case VSL_ERROR_NULL_PTR:
        case VSL_ERROR_CPU_NOT_SUPPORTED:
        case VSL_RNG_ERROR_INVALID_BRNG_INDEX:
        case VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED:
        case VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED:
        case VSL_RNG_ERROR_SKIPAHEADEX_UNSUPPORTED:
        case VSL_RNG_ERROR_BRNGS_INCOMPATIBLE:
        case VSL_RNG_ERROR_BAD_STREAM:
        case VSL_RNG_ERROR_BRNG_TABLE_FULL:
        case VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE:
        case VSL_RNG_ERROR_BAD_WORD_SIZE:
        case VSL_RNG_ERROR_BAD_NSEEDS:
        case VSL_RNG_ERROR_BAD_NBITS:
        case VSL_RNG_ERROR_QRNG_PERIOD_ELAPSED:
        case VSL_RNG_ERROR_LEAPFROG_NSTREAMS_TOO_BIG:
        case VSL_RNG_ERROR_BRNG_NOT_SUPPORTED:
        case VSL_RNG_ERROR_BAD_UPDATE:
        case VSL_RNG_ERROR_NO_NUMBERS:
        case VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM:
        case VSL_RNG_ERROR_NONDETERM_NOT_SUPPORTED:
        case VSL_RNG_ERROR_NONDETERM_NRETRIES_EXCEEDED:
        case VSL_RNG_ERROR_ARS5_NOT_SUPPORTED:
        case VSL_DISTR_MULTINOMIAL_BAD_PROBABILITY_ARRAY:
        case VSL_RNG_ERROR_FILE_CLOSE:
        case VSL_RNG_ERROR_FILE_OPEN:
        case VSL_RNG_ERROR_FILE_WRITE:
        case VSL_RNG_ERROR_FILE_READ:
        case VSL_RNG_ERROR_BAD_FILE_FORMAT:
        case VSL_RNG_ERROR_UNSUPPORTED_FILE_VER:
        case VSL_RNG_ERROR_BAD_MEM_FORMAT:
            return "(TODO)"; //TODO
        default:
            return "Bad error code";
        }
    }
}

#endif
