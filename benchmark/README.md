# Physica Benchmark

*Physica Benchmarks* are used for monitoring performance regression. Since performance testing is highly dependent on the platform, the test results are for reference only.$^{[1]}$。

## Notes

对于随机数质量要求不高的非热点函数，生成器统一采用线性同余生成器(MCG);  
使用[[clang::noinline]]避免感兴趣的callee内联, 便于后续LLVM IR层次分析;  

## References

[1] Physica Benchmark Results; <https://gitee.com/newsigma/llvm-opt-benchmark/tree/main/bench/Physica/Results>
