#############################################Settings#############################################
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
##############################################Libs################################################
find_package(Qt5 COMPONENTS Core Gui Charts REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_GUI *.cpp *.h)
add_library(PhysicaGui SHARED ${SOURCES_GUI})
#############################################Options##############################################
target_link_libraries(PhysicaGui PhysicaCore PhysicaInterpreter Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Charts)