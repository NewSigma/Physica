# ML
physica_add_example(DecisionTree ML/DecisionTree.cpp)
physica_add_example(Mnist ML/Mnist.cpp)
if(${PHYSICA_CUDA})
    physica_add_example(Mnist_cuda ML/Mnist.cu)
endif()

add_subdirectory(Math)
add_subdirectory(Physics)
