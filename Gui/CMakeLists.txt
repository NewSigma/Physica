#############################################Settings#############################################
set(CMAKE_PREFIX_PATH "~/Qt5.13.2/5.13.2/gcc_64/lib/cmake")
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
##############################################Libs################################################
find_package(Qt5 COMPONENTS Core Gui Widgets REQUIRED)
#############################################Project##############################################
if (TARGET Qt5::Core)
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_POSITION_INDEPENDENT_CODE "ON")
    set(CMAKE_CXX_COMPILE_OPTIONS_PIE "-fPIC")
endif()

file(GLOB_RECURSE SOURCES_GUI *.cpp *.h)
add_library(Physica_Gui ${SOURCES_GUI})
#############################################Options##############################################
target_link_libraries(Physica_Gui Qt5::Core Qt5::Gui Qt5::Widgets)