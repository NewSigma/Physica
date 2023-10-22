include(../cmake/Modules/PhysicaAddBenchmark.cmake)

# RPMD
physica_add_benchmark(HardCore RPMD/HardCore.cpp)
physica_add_benchmark(ParaH RPMD/ParaH.cpp)
if(${PHYSICA_CUDA})
    physica_add_benchmark(HardCore_cuda RPMD/HardCore.cu)
    physica_add_benchmark(ParaH_cuda RPMD/ParaH.cu)
    if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug") # nvcc 12.2 internal compiler error in debug mode
        physica_add_benchmark(ParaH_auto RPMD/ParaH_auto.cu)
    endif()
endif()
# Anonymous
physica_add_benchmark(KSSolver KSSolver.cpp)
physica_add_benchmark(Q_TIP4P Q_TIP4P.cpp)
physica_add_benchmark(Q_TIP4P2D Q_TIP4P2D.cpp)
