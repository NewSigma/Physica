# AI
physica_add_example(DecisionTree AI/DecisionTree.cpp)
physica_add_example(Mnist AI/Mnist.cpp)
if(${PHYSICA_CUDA})
    physica_add_example(Mnist_cuda AI/Mnist.cu)
endif()
# Math
physica_add_example(DuffingOscillator Math/DuffingOscillator.cpp)
physica_add_example(FokkerPlanck Math/FokkerPlanck.cpp)
physica_add_example(ProbDistribution2D Math/ProbDistribution2D.cpp)

add_subdirectory(Physics)
