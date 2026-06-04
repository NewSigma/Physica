add_library(PhysicaCore)

add_subdirectory(Exception)
add_subdirectory(IO)
add_subdirectory(Math)
add_subdirectory(Parallel)
add_subdirectory(Physics)
add_subdirectory(Scalar)
add_subdirectory(Utils)

target_sources(PhysicaCore PRIVATE Version.cpp)
target_sources(PhysicaCore PUBLIC Malloc.cpp) # Declared as PUBLIC, other targets will automatically overload operator new/delete
target_link_libraries(PhysicaCore PUBLIC Physica::Core)

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
install(FILES Malloc.cpp DESTINATION share)
