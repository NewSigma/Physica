# Physica

Physica is a high-performance and scalable template library for computational physics, dedicated to providing a rapid iteration platform for new algorithms. It is not a textbook or a ready-to-use solution. Physica provides several basic components and does not restrict creativity. It helps users implement ideas in less than 200 lines of code.

Physica is also an open source platform that maintains scientific code written in Physica and related data. We expect Physica, as a platform and not just a software, to promote the reproducibility of scientific results. The open source community continuously ensures that the results remain readable, reproducible and accurate.

## Table of Contents

- [Design philosophy](#philosophy)
- [Features](#feature)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Design philosophy <a id="philosophy"></a>

1. Be complex, strict and precise

2. Optimal performance and scalability: Balancing between performance, scalability and ease of use is a trade-off. Physica targets on performance and scalability, taking ease of use as a secondary concern if necessary.

## Features <a id="feature"></a>

- A linear algebra library that leverages SIMD and GPU acceleration, offers automatic differentiation support and multiprecision capabilities
- Basic math library: Provides commonly used functions such as ODE, PDE, special functions, optimization, statistics, etc.
- Classical molecular dynamics(MD), path integral molecular dynamics(PIMD) simulations using NVE, NVT and NPT ensembles
- Multithreads and CUDA parallel support
- 2D and 3D plotting support
- Template meta algorithm: optimal algorithm and parallelism strategy selection at compiling time
- Useful components that facilitate the construction of efficient computational physics programs

## Usage

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
[pybind11](https://github.com/pybind/pybind11) 2.11 (Optional, Python binding, Experimental)  
[LLVM](https://llvm.org/) 17.0.6 (Optional, Python binding, Experimental)  
[libffi](https://github.com/libffi/libffi/) 3.4.6 (Optional, Python binding, Experimental)  
[vectorclass](https://github.com/vectorclass/version2) 2.01.03 (Bundled)  

Compile Physica using the following command:  

```
mkdir -p /path/to/Physica/build
cd /path/to/Physica/build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install -DCMAKE_BUILD_TYPE=Release ../ # Note: Set CMake options as you need, refer to CMakeLists.txt for a full list of options
make install -j<N>
```

Test Physica using the following command:  

```
cd /path/to/Physica/build/test
ctest -j<N>
```

Do not use Physica if any test were failed.  

To use Physica, examples that under /path/to/install/Physica/examples are ready to use(Note: BUILD_EXAMPLES is OFF by default). Link Physica as a part of your project for production use.

## Maintainers

[@NewSigma](NewSigma@163.com)

## Contributing

Feel free to dive in! Open an issue or submit PRs.

### Something you can do

1.Test Physica on your machine and report bugs.  

2.Tell us new features you want.  

3.Fix bugs, add new features, tests and examples.  

4.Reproduce papers and contribute your results to the community  

5.Improve the documentation.  

### Contributors

## License

[GPLv3](LICENSE) © Weibo He
