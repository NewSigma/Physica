include(../cmake/Modules/PhysicaAddExample.cmake)

physica_add_example(CheckGPU CheckGPU.cu)

physica_add_example(CrossSection CrossSection/CrossSection.cpp)
physica_add_example(WaterDrop WaterDrop/WaterDrop.cpp)
physica_add_example(DensityPlotExample DensityPlotExample.cpp)
physica_add_example(DuffingOscillator DuffingOscillator.cpp)
physica_add_example(Histogram Histogram.cpp)
physica_add_example(StochasticRungeKutta StochasticRungeKutta.cpp)
