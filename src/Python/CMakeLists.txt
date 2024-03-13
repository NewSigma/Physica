##############################################Libs################################################
find_package(LLVM REQUIRED CONFIG)
find_package(Clang REQUIRED CONFIG)

if (NOT ${LLVM_ENABLE_RTTI})
    message(FATAL_ERROR "Build LLVM with RTTI is required")
endif()
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_PYTHON *.cpp)
if(${PHYSICA_CUDA})
    file(GLOB_RECURSE CUDA_SOURCES_PYTHON *.cu)
    list(APPEND SOURCES_PYTHON ${CUDA_SOURCES_PYTHON})
endif()

pybind11_add_module(PhysicaPython ${SOURCES_PYTHON})
#############################################Options##############################################
set(LLVM_INCLUDE_DIR "${LLVM_INSTALL_PREFIX}/lib/clang/${LLVM_VERSION_MAJOR}/include")
set(PHYSICA_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
target_include_directories(PhysicaPython SYSTEM PRIVATE ${LLVM_INCLUDE_DIRS})
target_compile_definitions(PhysicaPython PRIVATE ${LLVM_DEFINITIONS})
target_compile_definitions(PhysicaPython PRIVATE -DLLVM_RESOURCE_DIR="-resource-dir=${LLVM_LIBRARY_DIR}/clang/${LLVM_VERSION_MAJOR}")
target_compile_definitions(PhysicaPython PRIVATE -DPHYSICA_INCLUDE_DIR="-I${PHYSICA_INCLUDE_DIR}")
target_compile_definitions(PhysicaPython PRIVATE -DPHYSICA_CUDA_ARCHITECTURES="--offload-arch=sm_${PHYSICA_CUDA_ARCHITECTURES}")

target_link_libraries(PhysicaPython PRIVATE PhysicaCore PRIVATE clang-cpp)

install(TARGETS PhysicaPython DESTINATION lib)
