#############################################Settings#############################################
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
##############################################Libs################################################
find_package(Qt5 COMPONENTS Core Gui Widgets Charts REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_GUI *.*)
add_library(Physica_Gui SHARED ${SOURCES_GUI})
#############################################Options##############################################
target_link_libraries(Physica_Gui Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Charts)