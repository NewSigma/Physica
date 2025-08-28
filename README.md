# Physica

Physica is a high-performance and scalable C++ template library, dedicated to providing a rapid iteration platform for new algorithms. It provides the following two main features:

- Differentiable linear algebra library that leverages SIMD and GPU acceleration
- Domain specific composable modules for computational physics

Physica is also an open source platform that maintains scientific code written in Physica and related data. We expect Physica, as a platform and not just a software, to promote the development of open science. The open source community continuously ensures that the results remain readable, reproducible and reliable.

## Table of Contents

- [Design philosophy](#philosophy)
- [Features](#features)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Design philosophy

Optimal performance and scalability: Balancing between performance, scalability and ease of use is a trade-off. Physica targets on performance and scalability, taking ease of use as a secondary concern if necessary.

## Features

- Operator fusion using template expressions for both CPU and GPU  
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

which installs the required dependencies to build Physica. Refer to [Install.md](doc/Install.md) for instructions on how to build Physica.

## Maintainers

Weibo He (NewSigma@163.com)  

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
