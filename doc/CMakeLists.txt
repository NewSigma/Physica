# Doc-* target generates diagrams
function(physica_add_doc NAME SOURCES)
    add_executable(Doc-${NAME})
    target_sources(Doc-${NAME} PRIVATE ${SOURCES} ${ARGN})
    target_precompile_headers(Doc-${NAME} REUSE_FROM PhysicaCore)
    target_link_libraries(Doc-${NAME} PhysicaCore PhysicaGui ${GPerfTools_LIBRARY})
    set_target_properties(Doc-${NAME} PROPERTIES AUTOMOC ON)
endfunction()

add_subdirectory(Core)
