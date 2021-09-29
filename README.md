# Physica

Physica provides several basic modules in computational physics.

There will be a lot of work to do in the future

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Install

This project uses the following packages:  

[Qt](https://www.qt.io/)  (Optional)  
[CUDA](https://developer.nvidia.com/cuda-downloads)  (Optional)  

Then you can compile Physica using the following commands:

```sh
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/where/Physica/should/be/installed/to /path/to/Physica
$ make install    (alternatively $ make -j<N> install)
```

## Usage

You can run the executable file in your console:

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
