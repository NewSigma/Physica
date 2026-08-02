add_subdirectory(IO)
add_subdirectory(Math)
add_subdirectory(Scalar)
## ML
physica_add_test(AP ML/AP.cpp)
physica_add_test(BayesOpt ML/BayesOpt.cpp)
physica_add_test(GPModel ML/GPModel.cpp)

add_subdirectory(Parallel)
## Utils
physica_add_test(Array Utils/Array.cpp)
physica_add_test(Range Utils/Range.cpp)

add_subdirectory(Physics)
