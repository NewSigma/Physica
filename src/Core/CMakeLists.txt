##############################################Libs################################################
#############################################Project##############################################
add_library(PhysicaCore)

file(GLOB_RECURSE HEADERS_CORE ${CMAKE_SOURCE_DIR}/include/Physica/Core/*.h)
file(GLOB_RECURSE SOURCES_CORE *.cpp)
if (NOT ${PHYSICA_HDF5})
    file(GLOB_RECURSE HEADERS_HDF5 ${CMAKE_SOURCE_DIR}/include/Physica/Core/IO/HDF5/*.h)
    list(REMOVE_ITEM HEADERS_CORE ${HEADERS_HDF5})

    file(GLOB_RECURSE SOURCES_HDF5 IO/HDF5/*.cpp)
    list(REMOVE_ITEM SOURCES_CORE ${SOURCES_HDF5})
endif()

if(${PHYSICA_CUDA})
    file(GLOB_RECURSE SOURCES_CUDA *.cu)
    list(APPEND SOURCES_CORE ${SOURCES_CUDA})
endif()

foreach(header ${HEADERS_CORE})
    string(FIND ${header} "Impl" pos)
    if (NOT ${pos} MATCHES -1)
        list(REMOVE_ITEM HEADERS_CORE ${header})
    endif()
endforeach()

if(NOT ${CMAKE_DISABLE_PRECOMPILE_HEADERS})
    target_compile_definitions(PhysicaCore PUBLIC PhysicaCore_EXPORTS)
    target_precompile_headers(PhysicaCore PRIVATE ${HEADERS_CORE})
endif()
target_sources(PhysicaCore PRIVATE ${SOURCES_CORE})

target_link_libraries(PhysicaCore PUBLIC ${FFTW3_LIBRARIES} PUBLIC pthread)

if(${PHYSICA_HDF5})
    target_link_libraries(PhysicaCore PUBLIC ${HDF5_CXX_SHARED_LIBRARY})
endif()

if(${PHYSICA_CUDA})
    target_link_libraries(PhysicaCore PUBLIC CUDA::cublas)
endif()

if(${PHYSICA_MKL})
    target_link_libraries(PhysicaCore PUBLIC MKL::MKL)
endif()

if(${PHYSICA_MPI})
    target_link_libraries(PhysicaCore PUBLIC MPI::MPI_CXX)
endif()

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
