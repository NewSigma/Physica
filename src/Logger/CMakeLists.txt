# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2020-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_library(PhysicaLogger)

add_subdirectory(Logger)

target_sources(PhysicaLogger PRIVATE
    LogBuffer.cpp
    LoggerRuntime.cpp
    LoggerTimer.cpp
)

set_target_properties(PhysicaLogger PROPERTIES
    EXPORT_NAME Logger
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

target_link_libraries(PhysicaLogger PUBLIC PhysicaCore)

install(TARGETS PhysicaLogger EXPORT PhysicaLoggerTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${Physica_INCLUDE_DIR}/Physica/Logger DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Physica)

physica_add_header_check(PhysicaLogger)
