# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2025 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
find_path(FFI_INCLUDE_DIR ffi.h)
find_library(FFI_LIBRARY NAMES libffi_pic.a libffi.a)

set(FFI_INCLUDE_DIRS ${FFI_INCLUDE_DIR})
set(FFI_LIBRARIES ${FFI_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFI DEFAULT_MSG FFI_INCLUDE_DIR FFI_LIBRARY)

mark_as_advanced(FFI_INCLUDE_DIR FFI_LIBRARY)
