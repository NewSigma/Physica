find_path(FFI_INCLUDE_DIR ffi.h)
find_path(FFI_LIBRARY_DIR libffi.a PATH_SUFFIXES lib)
find_library(FFI_LIBRARY libffi.a)

set(FFI_INCLUDE_DIRS ${FFI_INCLUDE_DIR})
set(FFI_LIBRARY_DIRS ${FFI_LIBRARY_DIR})
set(FFI_LIBRARIES ${FFI_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFI DEFAULT_MSG FFI_LIBRARY FFI_INCLUDE_DIR FFI_LIBRARY_DIR)

mark_as_advanced(FFI_INCLUDE_DIR FFI_LIBRARY_DIR FFI_LIBRARY)
