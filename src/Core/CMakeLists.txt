##############################################Libs################################################
find_package(Boost REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.cpp *.h)
add_library(PhysicaCore SHARED ${SOURCES_CORE})
#############################################Options##############################################
target_link_libraries(PhysicaCore Qt5::Core PhysicaUtils)
#Qt::Core will pass wrong option -fPIC to nvcc, so we have to remove it.
if(TARGET Qt5::Core)
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
endif()