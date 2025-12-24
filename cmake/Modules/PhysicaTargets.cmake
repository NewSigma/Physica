##############################################Core################################################
add_library(Physica::Core INTERFACE IMPORTED)

target_link_libraries(Physica::Core INTERFACE ${PhysicaCore_LIBRARIES})
target_link_libraries(Physica::Core INTERFACE ${FFTW3_LIBRARIES} INTERFACE pthread)

if(CMAKE_CXX_COMPILER_ID MATCHES IntelLLVM)
    target_link_libraries(Physica::Core INTERFACE IntelSYCL::SYCL_CXX)
    target_link_libraries(Physica::Core INTERFACE imf INTERFACE svml INTERFACE irng INTERFACE intlc)
endif()

if(${PHYSICA_HDF5})
    target_link_libraries(Physica::Core INTERFACE hdf5::hdf5 hdf5::hdf5_cpp)
endif()

if(${PHYSICA_CUDA})
    target_link_libraries(Physica::Core INTERFACE CUDA::cublas INTERFACE CUDA::curand INTERFACE CUDA::cusolver INTERFACE cudss)
endif()

if(${PHYSICA_MKL})
    target_link_libraries(Physica::Core INTERFACE MKL::MKL)
endif()

if(${PHYSICA_MPI})
    target_link_libraries(Physica::Core INTERFACE MPI::MPI_CXX)
endif()

if(${PHYSICA_MIMALLOC})
    target_link_libraries(Physica::Core INTERFACE mimalloc)
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
