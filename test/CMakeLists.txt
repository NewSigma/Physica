include(../cmake/Modules/PhysicaAddTest.cmake)
if(${BUILD_TESTING}) #Declared by CTest
    message(STATUS "Configuring tests.")
    physica_add_test(ConstTest ConstTest.cpp)
else()
    message(STATUS "Test disabled.")
endif()