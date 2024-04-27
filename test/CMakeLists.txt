include(CTest)
include(../cmake/Modules/PhysicaAddTest.cmake)

message(STATUS "Configuring tests.")

add_subdirectory(Core)
physica_add_test(Logger Logger/LoggerTest.cpp)
physica_add_test(DirStack Utils/DirStackTest.cpp)
