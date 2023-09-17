# Physica

Physica is a library that provides several basic modules in computational physics.

## Table of Contents

- [Design philosophy](#philosophy)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Design philosophy

1. Optimal performance and scalability: Balancing between performance, scalability and ease of use is a trade-off. In Physica, we target on performance and scalability. Ease of use is taken secondary.

2. Zero overhead abstraction: Users shall not pay for what they do not need. Physica is mainly composed of several sets of header files. Codes that users do not need will not be compiled.

3. Template meta algorithm: No seperate input file. Simulation parameters are declared as compile time constant. The program will choose best algorithm at compile time automatically making use of powerful C++ templates.

## Usage

Requirements:  

[fftw](http://www.fftw.org)  3.3.10  
[HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.14.1 (Optional)  
[Qt](https://www.qt.io/)  6.2.1  (Optional)  
[oneMKL](https://www.intel.com/) 2023.2.0 (Optional)
[CUDA](https://developer.nvidia.com/cuda-downloads)  12.0  (Optional)

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
