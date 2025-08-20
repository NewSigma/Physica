add_library(PhysicaCore)

add_subdirectory(Exception)
add_subdirectory(IO)

file(GLOB_RECURSE SOURCES_CORE Math/*.cpp Parallel/*.cpp Physics/*.cpp Scalar/*.cpp Utils/*.cpp)
target_sources(PhysicaCore PRIVATE ${SOURCES_CORE})

if(${PHYSICA_CUDA})
    file(GLOB_RECURSE SOURCES_CUDA Math/*.cu Parallel/*.cu Physics/*.cu Scalar/*.cu Utils/*.cu)
    target_sources(PhysicaCore PRIVATE ${SOURCES_CUDA})
endif()

target_sources(PhysicaCore PRIVATE Version.cpp)
target_sources(PhysicaCore PUBLIC Malloc.cpp) # Declared as PUBLIC, other targets will automatically overload operator new/delete
target_link_libraries(PhysicaCore PUBLIC Physica::Core)

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
