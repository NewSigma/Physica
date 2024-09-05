find_path(FFTW3_INCLUDE_DIRS fftw3.h HINTS ${FFTW3_ROOT})
find_library(FFTW3_LIBRARY libfftw3.so HINTS ${FFTW3_ROOT})
find_library(FFTW3f_LIBRARY libfftw3f.so HINTS ${FFTW3_ROOT})

set(FFTW3_LIBRARIES ${FFTW3_LIBRARY} ${FFTW3f_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)

mark_as_advanced(FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)
