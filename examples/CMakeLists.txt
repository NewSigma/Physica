include(../cmake/Modules/PhysicaAddExample.cmake)

add_subdirectory(Core)
add_subdirectory(Gui)
# Annoymous
physica_add_example(CrossSection CrossSection/CrossSection.cpp)
physica_add_example(WaterDrop WaterDrop/WaterDrop.cpp)
physica_add_example(Ising Ising.cpp)
physica_add_example(VariationalSchrodinger VariationalSchrodinger.cpp)
