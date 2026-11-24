# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2020-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
add_library(PhysicaCore)

add_subdirectory(Exception)
add_subdirectory(IO)
add_subdirectory(Math)
add_subdirectory(Parallel)
add_subdirectory(Physics)
add_subdirectory(Scalar)
add_subdirectory(Utils)

include(GitHash.cmake)
target_sources(PhysicaCore PRIVATE Version.cpp ${CMAKE_CURRENT_BINARY_DIR}/GitHash.cpp)

set_target_properties(PhysicaCore PROPERTIES
    EXPORT_NAME Core
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

if(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_definitions(PhysicaCore PUBLIC _USE_MATH_DEFINES NOMINMAX)
endif()

target_include_directories(PhysicaCore PUBLIC $<BUILD_INTERFACE:${Physica_INCLUDE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(PhysicaCore PUBLIC pthread)
target_link_libraries(PhysicaCore PRIVATE FFTW3::fftw3 FFTW3::fftw3f)

if(CMAKE_CXX_COMPILER_ID MATCHES IntelLLVM)
    target_link_libraries(PhysicaCore PUBLIC IntelSYCL::SYCL_CXX)
    target_link_libraries(PhysicaCore PUBLIC imf svml irng intlc)
endif()

if(PHYSICA_PROFILE)
    target_include_directories(PhysicaCore SYSTEM PUBLIC $<BUILD_INTERFACE:${GPerfTools_INCLUDE_DIR}>)
endif()

if(PHYSICA_HDF5)
    target_compile_definitions(PhysicaCore PUBLIC PHYSICA_HDF5 H5_NO_DEPRECATED_SYMBOLS)
    target_link_libraries(PhysicaCore PUBLIC hdf5::hdf5)
endif()

if(PHYSICA_CUDA)
    target_compile_definitions(PhysicaCore PUBLIC PHYSICA_CUDA)
    target_link_libraries(PhysicaCore PUBLIC CUDA::cudart CUDA::cublas CUDA::curand CUDA::cusolver CUDA::nvrtc cudss)
endif()

if(PHYSICA_MKL)
    target_compile_definitions(PhysicaCore PUBLIC PHYSICA_MKL)
    target_link_libraries(PhysicaCore PUBLIC MKL::MKL)
endif()

if(PHYSICA_MPI)
    # Ubuntu 24.04 builds OpenMPI with CXX support. Drop it to avoid ambiguous name.
    target_compile_definitions(PhysicaCore PUBLIC PHYSICA_MPI OMPI_SKIP_MPICXX)
    target_link_libraries(PhysicaCore PUBLIC MPI::MPI_CXX)
endif()

if(PHYSICA_MIMALLOC)
    target_compile_definitions(PhysicaCore PUBLIC PHYSICA_MIMALLOC)
    target_link_libraries(PhysicaCore PUBLIC mimalloc)
    # Declared as PUBLIC, other targets will automatically overload operator new/delete
    target_sources(PhysicaCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Malloc.cpp> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DATADIR}/Physica/Malloc.cpp>)
    install(FILES Malloc.cpp DESTINATION ${CMAKE_INSTALL_DATADIR}/Physica)
endif()

if(PHYSICA_TRANSFORMS)
    target_compile_options(PhysicaCore INTERFACE $<BUILD_INTERFACE:-fpass-plugin=$<TARGET_FILE:PhysicaTransforms>>)
endif()

install(TARGETS PhysicaCore EXPORT PhysicaCoreTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${Physica_INCLUDE_DIR}/Physica/Core DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Physica)

physica_add_header_check(PhysicaCore)
