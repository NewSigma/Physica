# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2020-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_library(PhysicaGui)

add_subdirectory(Plot)

set_target_properties(PhysicaGui PROPERTIES
    EXPORT_NAME Gui
    INSTALL_RPATH \${ORIGIN}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

target_compile_definitions(PhysicaGui PUBLIC PHYSICA_GUI)
target_link_libraries(PhysicaGui PUBLIC PhysicaCore Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::DataVisualization Qt6::Svg)

install(TARGETS PhysicaGui EXPORT PhysicaGuiTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${Physica_INCLUDE_DIR}/Physica/Gui DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Physica)

physica_add_header_check(PhysicaGui)
