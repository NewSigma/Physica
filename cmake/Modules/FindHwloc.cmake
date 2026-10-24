find_path(HWLOC_INCLUDE_DIRS hwloc.h)
find_library(HWLOC_LIBRARY hwloc)

set(HWLOC_LIBRARIES ${HWLOC_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hwloc DEFAULT_MSG HWLOC_INCLUDE_DIRS HWLOC_LIBRARIES)

mark_as_advanced(HWLOC_INCLUDE_DIRS HWLOC_LIBRARY HWLOC_LIBRARIES)
