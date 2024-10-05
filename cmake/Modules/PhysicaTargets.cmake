add_library(Physica::Core INTERFACE IMPORTED)

target_include_directories(Physica::Core INTERFACE ${Physica_INCLUDE_DIR})
target_link_libraries(Physica::Core INTERFACE ${PhysicaCore_LIBRARIES})
