##############################################Libs################################################
#############################################Project##############################################
add_library(PhysicaCore)

file(GLOB_RECURSE SOURCES_CORE *.cpp)
if (NOT ${PHYSICA_HDF5})
    file(GLOB_RECURSE SOURCES_HDF5 IO/HDF5/*.cpp)
    list(REMOVE_ITEM SOURCES_CORE ${SOURCES_HDF5})
endif()

if(${PHYSICA_CUDA})
    file(GLOB_RECURSE CUDA_SOURCES_CORE *.cu)
    list(APPEND SOURCES_CORE ${CUDA_SOURCES_CORE})
endif()
target_sources(PhysicaCore PRIVATE ${SOURCES_CORE})

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaCore ${FFTW3_LIBRARIES} ${HDF5_CXX_SHARED_LIBRARY} pthread)

if(${PHYSICA_CUDA})
    target_link_libraries(PhysicaCore CUDA::cublas)
endif()

if(${PHYSICA_MKL})
    target_link_libraries(PhysicaCore MKL::MKL)
endif()

if(${PHYSICA_MPI})
    target_link_libraries(PhysicaCore MPI::MPI_CXX)
endif()
