#############################################Settings#############################################
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)

if(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    add_compile_options(/Wall)
    add_compile_options(/wd4068 /wd4201 /wd4242 /wd4244 /wd4251 /wd4267 /wd4275 /wd4305 /wd4365 /wd4458 /wd4514 /wd4623 /wd4625 /wd4626 /wd4710 /wd4711 /wd4800 /wd4819 /wd4820 /wd4996 /wd5026 /wd5027 /wd5045 /wd5051 /wd5219)
    add_definitions(-D_USE_MATH_DEFINES)
    add_definitions(-DNOMINMAX)
else()
    add_compile_options(-Wall -Wfatal-errors -mrdrnd -march=native -fno-semantic-interposition -fno-plt -fno-math-errno -fno-trapping-math -fno-signed-zeros)
    add_link_options(-Wl,-Bsymbolic,--as-needed,--gc-sections)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Og>)
    elseif(CMAKE_BUILD_TYPE MATCHES Release)
        add_link_options(-Wl,--strip-all)
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
        add_compile_options(-Wextra)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang OR CMAKE_CXX_COMPILER_ID MATCHES IntelLLVM)
        add_compile_options(-fassume-sane-operator-new -fassume-nothrow-exception-dtor)
        # Workaround for P2014R0
        #
        # Reference:
        # [1] https://github.com/llvm/llvm-project/issues/56671
        add_compile_options(-fcoro-aligned-allocation)
        add_link_options(-lstdc++) # Add this if you prefer libstdc++

        if (CMAKE_CXX_COMPILER_ID MATCHES IntelLLVM)
            find_package(IntelSYCL REQUIRED)
            get_target_property(IntelSYCL_LIBRARY_PATH IntelSYCL::SYCL_CXX INTERFACE_LINK_DIRECTORIES)
            set(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} ${IntelSYCL_LIBRARY_PATH})
            set(CMAKE_INSTALL_RPATH ${CMAKE_BUILD_RPATH} ${IntelSYCL_LIBRARY_PATH})
        endif()
    else()
        message(FATAL_ERROR "Unknown compiler")
    endif()

    if(FALSE)
        add_compile_options(-fsanitize=address)
        add_link_options(-fsanitize=address)
    endif()
endif()

if(${PHYSICA_CUDA})
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
    set(CMAKE_CUDA_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED})
    set(CMAKE_CUDA_VISIBILITY_PRESET ${CMAKE_CXX_VISIBILITY_PRESET})
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_definitions(-DPHYSICA_CUDA)

    if(CMAKE_CUDA_COMPILER_ID MATCHES Clang)
        add_compile_options(-Wno-unknown-cuda-version)
    else()
        # clangd does not work with response file
        # Reference: https://github.com/clangd/clangd/discussions/1676
        set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_INCLUDES 0)

        set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -Xcompiler -Og")
        set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO} -G -dopt=on")
        # Warning 20011: Call host function from host-device function
        # Warning 20208: Use long double in device code
        set(CMAKE_CUDA_FLAGS 
            ${CMAKE_CUDA_FLAGS}
            --expt-relaxed-constexpr
            --extended-lambda
            --device-entity-has-hidden-visibility true
            --ftz true
            --default-stream per-thread
            --Wreorder
            --Wdefault-stream-launch
            --Wext-lambda-captures-this
            --Wno-deprecated-gpu-targets
            --diag-suppress 20011
            --diag-suppress 20208
            ${CMAKE_CXX_FLAGS})
        string(REPLACE ";" " " CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
    endif()

    set(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} ${CUDAToolkit_LIBRARY_DIR})
    set(CMAKE_INSTALL_RPATH ${CMAKE_BUILD_RPATH} ${CUDAToolkit_LIBRARY_DIR})
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT Result OUTPUT Output)
    if(${Result})
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(WARNING "LTO Check failed, ignoring")
        message(TRACE ${Output})
    endif()
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
##############################################Libs################################################
# FFTW3
find_package(FFTW3 REQUIRED)
include_directories(SYSTEM ${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})

if (${PHYSICA_HDF5})
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    add_definitions(-DPHYSICA_HDF5 -DH5_NO_DEPRECATED_SYMBOLS)
endif()

if(${PHYSICA_PROFILE})
    find_package(GPerfTools REQUIRED)
    include_directories(SYSTEM ${GPerfTools_INCLUDE_DIR})
endif()

if(${PHYSICA_GUI})
    find_package(Qt6 COMPONENTS Core Gui Widgets Charts DataVisualization Svg REQUIRED)
    add_definitions(-DPHYSICA_GUI)
endif()

if(${PHYSICA_MKL})
    set(MKL_ARCH intel64)
    set(MKL_LINK dynamic)
    set(MKL_THREADING sequential)
    set(MKL_INTERFACE_FULL intel_ilp64)
    find_package(MKL REQUIRED)
    add_definitions(-DPHYSICA_MKL)
    include_directories(SYSTEM ${MKL_INCLUDE})
endif()

if(${PHYSICA_MPI})
    set(MPI_CXX_COMPILER_FLAGS -cc=${CMAKE_CXX_COMPILER})
    set(MPI_CXX_COMPILE_OPTIONS "")
    find_package(MPI REQUIRED)
    add_definitions(-DPHYSICA_MPI)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
endif()

if(${PHYSICA_CUDA})
    include_directories(SYSTEM ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    set(CUDA_MATH_INCLUDE_DIRECTORIES ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}/../../../../../math_libs/include)
    if (EXISTS ${CUDA_MATH_INCLUDE_DIRECTORIES}) # NVIDIA HPC SDK seems do not put math includes in standard location
        include_directories(SYSTEM ${CUDA_MATH_INCLUDE_DIRECTORIES})
    endif()
endif()
