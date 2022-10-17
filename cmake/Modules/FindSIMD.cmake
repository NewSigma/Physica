# Reference:
# [1] pytorch (https://github.com/pytorch/pytorch/blob/master/cmake/Modules/FindAVX.cmake)
include(CheckCSourceRuns)
include(CheckCXXSourceRuns)

SET(AVX_CODE "
    #include <immintrin.h>
    int main()
    {
        __m256 a;
        a = _mm256_set1_ps(0);
        return 0;
    }
")

SET(AVX2_CODE "
    #include <immintrin.h>
    int main()
    {
        __m256i a = {0};
        a = _mm256_abs_epi16(a);
        __m256i x;
        _mm256_extract_epi64(x, 0); // we rely on this in our AVX2 code
        return 0;
    }
")

SET(AVX512_CODE "
    #include <immintrin.h>
    int main()
    {
        __m512i a = _mm512_set_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0);
        __m512i b = a;
        __mmask64 equality_mask = _mm512_cmp_epi8_mask(a, b, _MM_CMPINT_EQ);
        return 0;
    }
")

macro(physica_check_simd type flags)
    set(__FLAG_I 1)
    set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
    foreach(__FLAG ${flags})
        if(NOT CXX_${type}_FOUND)
            set(CMAKE_REQUIRED_FLAGS ${__FLAG})
            CHECK_CXX_SOURCE_RUNS("${${type}_CODE}" CXX_HAS_${type}_${__FLAG_I})
            if(CXX_HAS_${type}_${__FLAG_I})
                set(CXX_${type}_FOUND TRUE CACHE BOOL "CXX ${type} support")
                set(CXX_${type}_FLAGS "${__FLAG}" CACHE STRING "CXX ${type} flags")
            endif()
            MATH(EXPR __FLAG_I "${__FLAG_I}+1")
        endif()
    endforeach()
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

    if(NOT CXX_${type}_FOUND)
        set(CXX_${type}_FOUND FALSE CACHE BOOL "CXX ${type} support")
        set(CXX_${type}_FLAGS "" CACHE STRING "CXX ${type} flags")
    endif()

    mark_as_advanced(CXX_${type}_FOUND CXX_${type}_FLAGS)
    add_definitions(-DPHYSICA_${type})
endmacro()

physica_check_simd("AVX" " ;-mavx;/arch:AVX")
physica_check_simd("AVX2" " ;-mavx2 -mfma;/arch:AVX2")
physica_check_simd("AVX512" " ;-mavx512f -mavx512dq -mavx512vl -mavx512bw -mfma;/arch:AVX512")
