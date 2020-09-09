include(../cmake/Modules/PhysicaAddTest.cmake)
if(${BUILD_TESTING}) #Declared by CTest
    message(STATUS "Configuring tests.")
    physica_add_test(ConstTest Core/MultiPrecision/ConstTest.cpp)
    physica_add_test(ScalarTest Core/MultiPrecision/ScalarTest.cpp)
else()
    message(STATUS "Test disabled.")
endif()