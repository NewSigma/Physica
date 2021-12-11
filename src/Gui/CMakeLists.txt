#############################################Settings#############################################
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
##############################################Libs################################################
find_package(Qt6 COMPONENTS Core Gui Widgets Charts DataVisualization REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_GUI *.cpp *.h)
add_library(PhysicaGui SHARED ${SOURCES_GUI})

install(TARGETS PhysicaGui DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Gui DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaGui PhysicaCore PhysicaInterpreter Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::DataVisualization)