add_library(PhysicaCore)

file(GLOB_RECURSE SOURCES_CORE *.cpp)
if (NOT ${PHYSICA_HDF5})
    file(GLOB_RECURSE SOURCES_HDF5 IO/HDF5/*.cpp)
    list(REMOVE_ITEM SOURCES_CORE ${SOURCES_HDF5})
endif()

if(${PHYSICA_CUDA})
    file(GLOB_RECURSE SOURCES_CUDA *.cu)
    list(APPEND SOURCES_CORE ${SOURCES_CUDA})
endif()

target_sources(PhysicaCore PRIVATE ${SOURCES_CORE})
target_link_libraries(PhysicaCore PUBLIC Physica::Core)

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
