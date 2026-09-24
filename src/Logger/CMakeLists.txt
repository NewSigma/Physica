add_library(PhysicaLogger)

add_subdirectory(Logger)

target_sources(PhysicaLogger PRIVATE
    LogBuffer.cpp
    LoggerRuntime.cpp
    LoggerTimer.cpp
)

target_link_libraries(PhysicaLogger PUBLIC Physica::Logger)

set_target_properties(PhysicaLogger PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

install(TARGETS PhysicaLogger DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Logger DESTINATION include/Physica)
