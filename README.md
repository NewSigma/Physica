# Physica

Physica is a high performance and scalable template library for computational physics which devotes to provide a rapid iteration platform for new algorithms. Physica is not a textbook nor ready to use solutions. It is expected users write their own codes and usually implement ideas in less than 200 lines of code.

## Table of Contents

- [Design philosophy](#philosophy)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Design philosophy

1. Optimal performance and scalability: Balancing between performance, scalability and ease of use is a trade-off. In Physica, we target on performance and scalability. The use of templates enables Physica face the challenge of heterogeneous computing. Ease of use is taken secondary if necessary.

2. Zero overhead abstraction: Users shall not pay for what they do not need. Physica is mainly composed of several sets of header files. Codes that users do not need will not be compiled.

3. Template meta algorithm: No seperate input file. Simulation parameters are declared as compile time constant. The program will choose best algorithm and parallization strategy at compile time automatically making use of powerful C++ templates.

4. Self explaination: As little documents and comments as necessary. Documents and comments often lag behind codes. The Physica API are designed to match users' intuition.

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
