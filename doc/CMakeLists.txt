# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2025 Weibo He
# SPDX-License-Identifier: MPL-2.0
# Doc-* target generates diagrams
function(physica_add_doc NAME SOURCES)
    add_executable(Doc-${NAME})
    target_sources(Doc-${NAME} PRIVATE ${SOURCES} ${ARGN})
    target_link_libraries(Doc-${NAME} PhysicaCore PhysicaGui ${GPerfTools_LIBRARY})
endfunction()

add_subdirectory(Core)
