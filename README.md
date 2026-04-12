# Physica

[![License](https://img.shields.io/badge/License-GPL_3.0-blue.svg)](https://opensource.org/license/gpl-3-0) [![Linux](https://img.shields.io/badge/os-Linux-green.svg)](https://shields.io/)

*Physica* is a high-performance and scalable C++ template library, dedicated to providing a rapid iteration platform for new algorithms. It provides the following two main features:

- Differentiable linear algebra library that leverages SIMD and GPU acceleration
- Domain specific composable modules for computational physics

It is also an open source platform that maintains scientific code and related data. We expect it as a platform and not just a software, to promote the development of open science. The open source community continuously ensures that the results remain readable, reproducible and reliable.

## Table of Contents

- [Features](#features)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Features

- Operator fusion using expression templates for both CPU and GPU  
- Coroutine based automatic differentiation implementation  
- Basic math library: Provides commonly used functions such as ODE, PDE, special functions, optimization, statistics, etc.
- Template meta algorithm: optimal algorithm and parallelism strategy selection at compiling time
- Multi-threads, MPI and CUDA parallel support
- 2D and 3D plotting support

## Usage

We provide a conda package to help with fast deployment:

``` Bash
conda install -c conda-forge -c nvidia newsigma::physica
```

which should install the required dependencies. Refer to [Install.md](doc/Install.md) for a detailed guide on using *Physica*.

## Maintainers

Weibo He (<NewSigma@163.com>)  

## Contributing

Feel free to dive in! Open an issue or submit PRs.

### Something you can do

1.Test it on your machine and report bugs.  

2.Tell us new features you want.  

3.Fix bugs, add new features, tests and examples.  

4.Reproduce papers and contribute your results to the community  

5.Improve the documentation.  

### Contributors

## License

[GPLv3](LICENSE) © Weibo He
