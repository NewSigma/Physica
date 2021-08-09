find_path(GPerfTools_INCLUDE_DIR gperftools/tcmalloc.h)

find_library(GPerfTools_LIBRARY libprofiler.so)

set(GPerfTools_LIBRARIES ${GPerfTools_LIBRARY})
set(GPerfTools_INCLUDE_DIRS ${GPerfTools_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPerfTools DEFAULT_MSG GPerfTools_LIBRARY GPerfTools_INCLUDE_DIR)

mark_as_advanced(GPerfTools_INCLUDE_DIR GPerfTools_LIBRARY)
