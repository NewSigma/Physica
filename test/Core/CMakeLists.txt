add_subdirectory(IO)
add_subdirectory(Math)
add_subdirectory(Scalar)
## ML
physica_add_test(AP ML/AP.cpp)

add_subdirectory(Parallel)
## Utils
physica_add_test(Array Utils/Array.cpp)

add_subdirectory(Physics)
