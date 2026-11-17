find_path(FFTW3_INCLUDE_DIRS fftw3.h HINTS ${FFTW3_ROOT})
find_library(FFTW3_LIBRARY fftw3 HINTS ${FFTW3_ROOT})
find_library(FFTW3f_LIBRARY fftw3f HINTS ${FFTW3_ROOT})

set(FFTW3_LIBRARIES ${FFTW3_LIBRARY} ${FFTW3f_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)

if(FFTW3_FOUND)
    if(NOT TARGET FFTW3::fftw3)
        add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::fftw3 PROPERTIES IMPORTED_LOCATION "${FFTW3_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()

    if(NOT TARGET FFTW3::fftw3f)
        add_library(FFTW3::fftw3f UNKNOWN IMPORTED)
        set_target_properties(FFTW3::fftw3f PROPERTIES IMPORTED_LOCATION "${FFTW3f_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()
endif()

mark_as_advanced(FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)
