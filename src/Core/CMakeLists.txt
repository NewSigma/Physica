add_library(PhysicaCore)

add_subdirectory(Exception)
add_subdirectory(IO)
add_subdirectory(Math)
add_subdirectory(Parallel)
add_subdirectory(Physics)
add_subdirectory(Scalar)
add_subdirectory(Utils)

target_sources(PhysicaCore PRIVATE Version.cpp)
target_link_libraries(PhysicaCore PUBLIC Physica::Core)

set_target_properties(PhysicaCore PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)

if(${PHYSICA_MIMALLOC})
    target_sources(PhysicaCore PUBLIC Malloc.cpp) # Declared as PUBLIC, other targets will automatically overload operator new/delete
    install(FILES Malloc.cpp DESTINATION share)
endif()
