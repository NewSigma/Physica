# Physica

Physica aims to be a computational physics engine. This project is created to offer a platform 

where students can improve their mathematics, physics and programming level.
  
There are a lot of work to do in the future.


## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Background

The final goal of Physica is to become a software like Matlab or Mathematica. This is a hard work

and we want to program to improve ourselves instead of finish the goal in the short term.

## Install

This project uses the following software. Go check them out if you don't have them locally installed:  

[CMake](https://cmake.org/)  
[Qt](https://www.qt.io/)  
[CUDA](https://developer.nvidia.com/cuda-downloads).

Then you can compile Physica using the following commands:

```sh
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/where/Physica/should/be/installed/to /path/to/Physica
$ make install    (alternatively $ make -j<N> install)
```

## Usage

You can run the executable file your console:

```sh
$ /path/where/Physica/should/be/installed/to/Physica
```

Or link Physica to your program.

## Maintainers

[@NewSigma](https://gitee.com/newsigma).

## Contributing

Feel free to dive in! Open an issue or submit PRs.

### Something you can do

1.Test Physica on your machine and report bugs.  

2.Tell us new features you want.  

3.Provide a better algorithm with mathematical prove.  

4.Fix bugs, add new features, tests and examples.  

5.Improve the documentation.  

### Contributors

## License

[GPLv3](LICENSE) © WeiBo He