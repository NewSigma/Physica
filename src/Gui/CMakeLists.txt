##############################################Libs################################################
find_package(Qt5 COMPONENTS Core Gui Widgets Charts REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_GUI *.cpp *.h)
add_library(Physica_Gui SHARED ${SOURCES_GUI})
#############################################Options##############################################
target_link_libraries(Physica_Gui Physica_Core Physica_Interpreter Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Charts)