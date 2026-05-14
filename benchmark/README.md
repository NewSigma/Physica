# Physica Benchmark

*Physica Benchmarks* is devoted to

- Select backend based on benchmark data
- Monitor performance regression

Since performance testing is highly dependent on the platform, the test results are for reference only.$^{[1]}$。

## Notes

For non-hotspot functions that do not require high-quality random numbers, the generator uniformly uses a multiplicative congruential generator (MCG).  
The [[clang::noinline]] attribute is used to prevent inlining of the callee of interest, facilitating subsequent LLVM IR-level analysis.  

## References

[1] Physica Benchmark Results; <https://gitee.com/newsigma/llvm-opt-benchmark/tree/main/bench/Physica/Results>
