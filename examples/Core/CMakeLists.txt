# AI
physica_add_example(DecisionTree AI/DecisionTree.cpp)
physica_add_example(Mnist AI/Mnist.cpp)
physica_add_example(Mnist_cuda AI/Mnist.cu)
# Math
physica_add_example(DuffingOscillator Math/DuffingOscillator.cpp)
physica_add_example(FokkerPlanck Math/FokkerPlanck.cpp)
physica_add_example(ProbabilityDistributionFunction2D Math/ProbabilityDistributionFunction2D.cpp)
# Physics
physica_add_example(HF Physics/ElectronicStructure/HF.cpp)
physica_add_example(ExactDiag Physics/ManyBody/ExactDiag.cpp)
physica_add_example(Correlation Physics/MD/Analyser/Correlation.cpp)
physica_add_example(RDF Physics/MD/Analyser/RDF.cpp)
physica_add_example(FireModel Physics/MD/KineticModel/FireModel.cpp)
physica_add_example(TodaModel Physics/MD/ForceModel/TodaModel.cpp)
if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug") # nvcc 12.2 internal compiler error in debug mode
    physica_add_example(HardCore_cuda Physics/MD/HardCore.cu)
endif()
physica_add_example(HardCore Physics/MD/HardCore.cpp)
physica_add_example(TestNVE Physics/MD/TestNVE.cpp)
