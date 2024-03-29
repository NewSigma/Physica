include(../cmake/Modules/PhysicaAddExample.cmake)

if(${PHYSICA_TORCH})
    physica_add_example(Regression AI/Regression.cpp)
endif()

if(${PHYSICA_CUDA})
    physica_add_example(CheckGPU CheckGPU.cu)
endif()

if(${PHYSICA_GUI})
    physica_add_example(CrossSection CrossSection/CrossSection.cpp)
    physica_add_example(WaterDrop WaterDrop/WaterDrop.cpp)
    physica_add_example(DensityPlotExample DensityPlotExample.cpp)
    physica_add_example(DuffingOscillator DuffingOscillator.cpp)
    physica_add_example(HF HF.cpp)
    physica_add_example(Histogram Histogram.cpp)
    physica_add_example(SurfacePlot SurfacePlot.cpp)
    physica_add_example(VariationalSchrodinger VariationalSchrodinger.cpp)
endif()

physica_add_example(Ising Ising.cpp)
