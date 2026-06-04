add_library(PhysicaLogger)

add_subdirectory(Logger)

target_sources(PhysicaLogger PRIVATE
    LogBuffer.cpp
    LoggerRuntime.cpp
    LoggerTimer.cpp
)

if(${PHYSICA_MIMALLOC})
    target_link_libraries(PhysicaLogger PRIVATE mimalloc)
endif()

install(TARGETS PhysicaLogger DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Logger DESTINATION include/Physica)
