# Benchmark

OS: Ubuntu24.04  
Compiler: clang 17.0.6 + NVCC 12.6.20  
Hardware: Intel(R) Core(TM) i7-14700KF + 32G + NVIDIA GeForce RTX 4060  
Date: 2024-11-17  
Maintainer: Weibo He (NewSigma@163.com)  

----------------------------------------------------------
Benchmark                Time             CPU   Iterations
----------------------------------------------------------
HubbardMatrix1D       29.9 s          29.9 s             1
TPQ                   1.14 s          1.14 s             1
HardCore             0.334 s         0.334 s             2
ParaH 108            0.164 s         0.089 s             7
ParaH 256            0.451 s         0.360 s             2
ParaH 500            0.363 s         0.298 s             2
ParaH 864            0.291 s         0.242 s             3
Q_TIP4P 1            0.685 s         0.685 s             1
Q_TIP4P 2             1.37 s          1.37 s             1
Q_TIP4P 4             1.37 s          1.37 s             1
HardCore_cuda         36.3 s          36.3 s             1
ParaH cuda 108       0.366 s         0.366 s             2
ParaH cuda 256       0.613 s         0.613 s             1
ParaH cuda 500       0.518 s         0.518 s             1
ParaH cuda 864       0.431 s         0.431 s             2
ParaH auto 108       0.200 s         0.108 s             5
ParaH auto 256       0.335 s         0.263 s             3
ParaH auto 500       0.273 s         0.239 s             3
ParaH auto 864       0.158 s         0.147 s             5
Q_TIP4P cuda         0.813 s         0.813 s             1
Mnist                 17.9 s          17.9 s             1
