enable_language(C)
##############################################Libs################################################
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)

find_package(nanobind REQUIRED CONFIG)
find_package(LLVM REQUIRED CONFIG)
find_package(Clang REQUIRED CONFIG)
find_package(FFI REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_PYTHON *.cpp)
if(${PHYSICA_CUDA})
    file(GLOB_RECURSE CUDA_SOURCES_PYTHON *.cu)
    list(APPEND SOURCES_PYTHON ${CUDA_SOURCES_PYTHON})
endif()
nanobind_add_module(PhysicaPython NOMINSIZE ${SOURCES_PYTHON})

set(LLVM_INCLUDE_DIR "${LLVM_INSTALL_PREFIX}/lib/clang/${LLVM_VERSION_MAJOR}/include")
target_include_directories(PhysicaPython SYSTEM PRIVATE ${LLVM_INCLUDE_DIRS})
target_include_directories(PhysicaPython SYSTEM PRIVATE ${FFI_INCLUDE_DIRS})

target_compile_definitions(PhysicaPython PRIVATE ${LLVM_DEFINITIONS})
target_compile_definitions(PhysicaPython PRIVATE -DLLVM_RESOURCE_DIR="${LLVM_LIBRARY_DIR}/clang/${LLVM_VERSION_MAJOR}")
target_compile_definitions(PhysicaPython PRIVATE -DHDF5_INCLUDE_DIR="${HDF5_INCLUDE_DIR_CPP}")
target_compile_definitions(PhysicaPython PRIVATE -DPHYSICA_CUDA_ARCHITECTURES="--offload-arch=sm_${PHYSICA_CUDA_ARCHITECTURES}")

target_link_libraries(PhysicaPython PRIVATE PhysicaCore PRIVATE clang-cpp PRIVATE ${FFI_LIBRARY})

set_target_properties(PhysicaPython PROPERTIES INSTALL_RPATH \${ORIGIN})

install(TARGETS PhysicaPython DESTINATION lib)
