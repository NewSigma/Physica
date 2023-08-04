include(../cmake/Modules/PhysicaAddBenchmark.cmake)

physica_add_benchmark(HardCore RPMD/HardCore.cpp)
physica_add_benchmark(HardCore_cuda RPMD/HardCore.cu)
physica_add_benchmark(KSSolver KSSolver.cpp)
physica_add_benchmark(ParaH ParaH.cpp)
physica_add_benchmark(Q_TIP4P Q_TIP4P.cpp)
physica_add_benchmark(Q_TIP4P2D Q_TIP4P2D.cpp)
