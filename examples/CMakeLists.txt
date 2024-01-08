include(../cmake/Modules/PhysicaAddExample.cmake)

if(${PHYSICA_GUI})
    # Core
    physica_add_example(DuffingOscillator Core/Math/DuffingOscillator.cpp)
    physica_add_example(ProbabilityDistributionFunction2D Core/Math/ProbabilityDistributionFunction2D.cpp)
    physica_add_example(HF Core/Physics/ElectronicStructure/HF.cpp)
    physica_add_example(RDF Core/Physics/MD/Analyser/RDF.cpp)
    physica_add_example(FireModel Core/Physics/MD/KineticModel/FireModel.cpp)
    physica_add_example(TodaModel Core/Physics/MD/ForceModel/TodaModel.cpp)
    if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug") # nvcc 12.2 internal compiler error in debug mode
        physica_add_example(HardCore_cuda Core/Physics/MD/HardCore.cu)
    endif()
    physica_add_example(HardCore Core/Physics/MD/HardCore.cpp)
    # Gui
    physica_add_example(DensityPlot Gui/DensityPlot.cpp)
    physica_add_example(GaussPlot Gui/GaussPlot.cpp)
    # Annoymous
    physica_add_example(CrossSection CrossSection/CrossSection.cpp)
    physica_add_example(WaterDrop WaterDrop/WaterDrop.cpp)
    physica_add_example(Ising Ising.cpp)
    physica_add_example(VariationalSchrodinger VariationalSchrodinger.cpp)
endif()

# Core
physica_add_example(DecisionTree Core/AI/DecisionTree.cpp)
physica_add_example(Correlation Core/Physics/MD/Correlation.cpp)
