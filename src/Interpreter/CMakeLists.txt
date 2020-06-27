##############################################Libs################################################
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_INTERPRETER *.cpp *.h)
add_library(Physica_Interpreter SHARED ${SOURCES_INTERPRETER})
#############################################Options##############################################
target_link_libraries(Physica_Interpreter Physica_Core Qt5::Core)