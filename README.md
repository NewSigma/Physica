# Physica

Physica is a high-performance and scalable template library for computational physics, dedicated to providing a rapid iteration platform for new algorithms. It is not a textbook or a ready-to-use solution. Rather, it is expected that users will write their own codes and implement ideas in less than 200 lines of code.

## Table of Contents

- [Design philosophy](#philosophy)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Design philosophy

1. Optimal performance and scalability: Balancing between performance, scalability and ease of use is a trade-off. Physica targets on performance and scalability, taking ease of use as a secondary concern if necessary.

2. Zero overhead abstraction: Users shall not pay for what they do not need. Physica is primarily composed of several sets of header files. Codes that are not required by users will not be compiled.

3. Template meta algorithm: No separate input file is required. Simulation parameters are declared as compile-time constants. The program automatically selects the best algorithm and parallelism strategy at compile time, leveraging powerful C++ templates.

4. Self explaination: Minimizing documentation and comments. Documentation and comments often lag behind code updates. The Physica API is designed to align with users' intuition.

## Usage

Requirements:  

[fftw](http://www.fftw.org)  3.3.10  
[HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.14.1  
[Qt](https://www.qt.io/)  6.2.1  (Optional, Plotting support)  
[CUDA](https://developer.nvidia.com/cuda-downloads)  12.0  (Optional)  
[oneMKL](https://www.intel.com/) 2023.2.0 (Optional)

To use Physica, simply compile and link Physica as a part of your project.

## Maintainers

[@NewSigma](https://gitee.com/newsigma).

## Contributing

Feel free to dive in! Open an issue or submit PRs.

### Something you can do

1.Test Physica on your machine and report bugs.  

2.Tell us new features you want.  

3.Provide a better implementation of algorithms.  

4.Fix bugs, add new features, tests and examples.  

5.Improve the documentation.  

### Contributors

## License

[GPLv3](LICENSE) © WeiBo He
