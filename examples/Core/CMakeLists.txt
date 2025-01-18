# AI
physica_add_example(DecisionTree AI/DecisionTree.cpp)
physica_add_example(Mnist AI/Mnist.cpp)
if(${PHYSICA_CUDA})
    #physica_add_example(Mnist_cuda AI/Mnist.cu) # NVCC 12.6 ICE
endif()

add_subdirectory(Math)
add_subdirectory(Physics)
