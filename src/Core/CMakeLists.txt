##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.cpp)
if(${PHYSICA_CUDA})
    file(GLOB_RECURSE CUDA_SOURCES_CORE *.cu)
    list(APPEND SOURCES_CORE ${CUDA_SOURCES_CORE})
endif()
add_library(PhysicaCore SHARED ${SOURCES_CORE})

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaCore PhysicaLogger PhysicaUtils ${FFTW3_LIBRARIES} ${FFTW3f_LIBRARIES} pthread)

if(${PHYSICA_MKL})
    target_link_libraries(PhysicaCore MKL::MKL)
endif()

if(${PHYSICA_HDF5})
    target_link_libraries(PhysicaCore ${HDF5_CXX_SHARED_LIBRARY})
endif()
