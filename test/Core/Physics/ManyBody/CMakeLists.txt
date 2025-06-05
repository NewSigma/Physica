physica_add_test(HubbardMatrix HubbardMatrix.cpp)
physica_add_test(KFermiRepr KFermiRepr.cpp)
physica_add_test(SpinState SpinState.cpp)
physica_add_test(TransIsingMatrix TransIsingMatrix.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    physica_add_test(DQMC DQMC.cpp)
endif()
