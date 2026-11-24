# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
find_package(Git QUIET REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_HASH
                ERROR_VARIABLE GIT_HASH_ERROR
                RESULT_VARIABLE GIT_HASH_RESULT
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT GIT_HASH_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to read the git revision (git returned ${GIT_HASH_RESULT}): ${GIT_HASH_ERROR}")
endif()
if (GIT_HASH STREQUAL "")
    message(FATAL_ERROR "Git revision is empty?")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GitHash.cpp.in GitHash.cpp @ONLY)
