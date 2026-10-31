# Physica Benchmark Agent Guide

## Conventions

1. Wrap the callee of interest with `PHYSICA_BENCH`; it expands to `[[clang::noinline]]` to keep the callee out of the benchmark, facilitating subsequent analysis at the LLVM IR level.
2. For non-hotspot functions that do not require high-quality random numbers, the generator uniformly uses a multiplicative congruential generator (MCG). Declare it as `using RandomSource = Random<MCG, std::mt19937::default_seed>;` and never spell the type out at the use site.
3. Use a fixed random seed and reset it before each benchmark run, since performance may be data-sensitive.
