##############################################Core################################################
add_library(Physica::Core INTERFACE IMPORTED)

target_link_libraries(Physica::Core INTERFACE ${PhysicaCore_LIBRARIES})
target_link_libraries(Physica::Core INTERFACE ${FFTW3_LIBRARIES} INTERFACE pthread)
target_include_directories(Physica::Core SYSTEM INTERFACE ${FFTW3_INCLUDE_DIRS})

if(CMAKE_CXX_COMPILER_ID MATCHES IntelLLVM)
    target_link_libraries(Physica::Core INTERFACE IntelSYCL::SYCL_CXX)
    target_link_libraries(Physica::Core INTERFACE imf INTERFACE svml INTERFACE irng INTERFACE intlc)
endif()

if(${PHYSICA_PROFILE})
    target_include_directories(Physica::Core SYSTEM INTERFACE ${GPerfTools_INCLUDE_DIR})
endif()

if(${PHYSICA_HDF5})
    target_link_libraries(Physica::Core INTERFACE hdf5::hdf5 hdf5::hdf5_cpp)
endif()

if(${PHYSICA_CUDA})
    target_link_libraries(Physica::Core INTERFACE CUDA::cublas INTERFACE CUDA::curand INTERFACE CUDA::cusolver INTERFACE cudss)
endif()

if(${PHYSICA_MKL})
    target_link_libraries(Physica::Core INTERFACE MKL::MKL)
    target_include_directories(Physica::Core SYSTEM INTERFACE ${MKL_INCLUDE})
endif()

if(${PHYSICA_MPI})
    target_link_libraries(Physica::Core INTERFACE MPI::MPI_CXX)
    target_include_directories(Physica::Core SYSTEM INTERFACE ${MPI_CXX_INCLUDE_DIRS})
endif()

if(${PHYSICA_MIMALLOC})
    target_link_libraries(Physica::Core INTERFACE mimalloc)
endif()

if(${PHYSICA_TRANSFORMS})
    target_compile_options(Physica::Core INTERFACE -fpass-plugin=$<TARGET_FILE:PhysicaTransforms>)
endif()
##############################################Gui################################################
if(${PHYSICA_GUI})
    add_library(Physica::Gui INTERFACE IMPORTED)

    target_link_libraries(Physica::Gui INTERFACE ${PhysicaGui_LIBRARIES})
    target_link_libraries(Physica::Gui INTERFACE PhysicaCore INTERFACE Qt6::Core INTERFACE Qt6::Gui INTERFACE Qt6::Widgets)
    target_link_libraries(Physica::Gui INTERFACE Qt6::Charts INTERFACE Qt6::DataVisualization INTERFACE Qt6::Svg)
endif()
##############################################Logger################################################
add_library(Physica::Logger INTERFACE IMPORTED)

target_link_libraries(Physica::Logger INTERFACE ${PhysicaLogger_LIBRARIES})
target_link_libraries(Physica::Logger INTERFACE PhysicaCore)
