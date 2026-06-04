add_subdirectory(PlotImpl)

target_sources(PhysicaGui PRIVATE
    DensityPlot.cpp
    GaussPlot.cpp
    MultiPlot.cpp
    Plot3D.cpp
    Plot.cpp
    PolarPlot.cpp
)
