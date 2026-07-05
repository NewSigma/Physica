physica_add_test(Gnuplot Gnuplot.cpp)
if(${PHYSICA_HDF5})
    physica_add_test(HDF5 HDF5.cpp)
endif()
