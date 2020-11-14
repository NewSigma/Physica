##############################################Libs################################################
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.cpp)
add_library(PhysicaUtils SHARED ${SOURCES_CORE})
#############################################Options##############################################
target_link_libraries(PhysicaUtils Qt5::Core)
