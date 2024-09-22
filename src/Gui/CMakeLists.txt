##############################################Libs################################################
find_package(Qt6 COMPONENTS Core Gui Widgets Charts DataVisualization REQUIRED)
#############################################Project##############################################
add_library(PhysicaGui)

set_target_properties(PhysicaGui PROPERTIES AUTOMOC ON)
set_target_properties(PhysicaGui PROPERTIES AUTOUIC ON)
set_target_properties(PhysicaGui PROPERTIES AUTORCC ON)

file(GLOB_RECURSE SOURCES_GUI *.cpp *.h)
target_sources(PhysicaGui PRIVATE ${SOURCES_GUI})
target_precompile_headers(PhysicaGui REUSE_FROM PhysicaCore)
target_link_libraries(PhysicaGui PhysicaCore Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::DataVisualization)

install(TARGETS PhysicaGui DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Gui DESTINATION include/Physica)
