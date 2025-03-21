# Install

We adopt CMake as our building system:

[CMake](https://cmake.org/) Refer CMakeLists.txt for minimum version  

A C++ compiler that support C++ 20, the following compilers passed our test:

[GCC](https://gcc.gnu.org/) N/A (Bug 104177)  
[clang](https://clang.llvm.org/) 17.0.6  
[IntelLLVM](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/dpc-compiler.html) 2025.0  

Addtional libraries:

[fftw](http://www.fftw.org)  3.3.10  
[HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.14.1  (Optional, Recommended data format)  
[oneMKL](https://www.intel.com/) 2024.2 (Optional)  
[CUDA](https://developer.nvidia.com/cuda-downloads)  12.8  (Optional)  
[Qt](https://www.qt.io/)  6.5.3  (Optional, Plotting support)  
[nanobind](https://github.com/wjakob/nanobind) 2.5 (Optional, Python binding, Concept Validation)  
[LLVM](https://llvm.org/) 17.0.6 (Optional, Python binding, Concept Validation)  
[libffi](https://github.com/libffi/libffi/) 3.4.6 (Optional, Python binding, Concept Validation)  
[vectorclass](https://github.com/vectorclass/version2) 2.01.03 (Bundled)  

Compile Physica using the following command:  

``` Bash
mkdir -p /path/to/Physica/build
cd /path/to/Physica/build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install -DCMAKE_BUILD_TYPE=Release ../ # Note: Set CMake options as you need, refer to CMakeLists.txt for a full list of options
cmake --build .
cmake --install .
```

Test Physica using the following command:  

``` Bash
cd /path/to/Physica/build/test
ctest -j<N>
```

Do not use Physica if any test were failed.  

To use Physica, examples that under /path/to/install/Physica/examples are ready to use(Note: BUILD_EXAMPLES is OFF by default). Link Physica as a part of your project for production use.
