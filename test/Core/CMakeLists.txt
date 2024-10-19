## AI
physica_add_test(AP AI/AP.cpp)
physica_add_test(Loss AI/Loss.cpp)

add_subdirectory(IO)
add_subdirectory(Math)
## MultiPrecision
physica_add_test(Const MultiPrecision/ConstTest.cpp)
physica_add_test(Diff MultiPrecision/Diff.cpp)
physica_add_test(Integer MultiPrecision/IntegerTest.cpp)
physica_add_test(Scalar MultiPrecision/ScalarTest.cpp)
## Parallel
physica_add_test(ThreadPool Parallel/ThreadPool.cpp)
## Utils
physica_add_test(DirStack Utils/DirStackTest.cpp)

add_subdirectory(Physics)
