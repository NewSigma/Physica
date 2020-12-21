##############################################Libs################################################
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_AI *.cpp *.h)
add_library(PhysicaAI SHARED ${SOURCES_AI})

install(TARGETS PhysicaAI DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/AI DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaAI PhysicaCore Qt5::Core)
#Qt::Core will pass wrong option -fPIC to nvcc, so we have to remove it.
if(TARGET Qt5::Core)
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
endif()