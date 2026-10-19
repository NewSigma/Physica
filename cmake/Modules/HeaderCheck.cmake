function(physica_collect_headers OUT_VAR)
    file(GLOB_RECURSE headers CONFIGURE_DEPENDS ${Physica_INCLUDE_DIR}/Physica/*.h)
    list(FILTER headers EXCLUDE REGEX "/[^/]*Impl[^/]*/")
    list(FILTER headers EXCLUDE REGEX "Impl\\.h$")
    # Headers depending on an optional dependency are out of reach unless it is enabled
    if(NOT PHYSICA_GUI)
        list(FILTER headers EXCLUDE REGEX "/Gui/")
    endif()
    if(NOT PHYSICA_HDF5)
        list(FILTER headers EXCLUDE REGEX "/IO/HDF5/|DeepModelDataset\\.h$")
    endif()
    if(NOT PHYSICA_MKL)
        list(FILTER headers EXCLUDE REGEX "_MKL\\.h$")
    endif()
    set(${OUT_VAR} ${headers} PARENT_SCOPE)
endfunction()

function(physica_add_header_check)
    foreach(dep IN LISTS ARGN)
        target_link_libraries(HeaderCheck INTERFACE ${dep})
    endforeach()
endfunction()

add_library(HeaderCheck INTERFACE)

physica_collect_headers(headers)
target_sources(HeaderCheck INTERFACE FILE_SET HEADERS BASE_DIRS ${Physica_INCLUDE_DIR} FILES ${headers})

set_target_properties(HeaderCheck PROPERTIES VERIFY_INTERFACE_HEADER_SETS TRUE)
