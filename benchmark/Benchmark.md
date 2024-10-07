# Benchmark

OS: Ubuntu24.04  
Compiler: GCC 9.5.0 + NVCC 12.5  
Hardware: Intel(R) Core(TM) i7-14700KF + 32G + NVIDIA GeForce RTX 4060  
Date: 2024-10-7  
Maintainer: Weibo He (NewSigma@163.com)  

----------------------------------------------------------
Benchmark                Time             CPU   Iterations
----------------------------------------------------------
HubbardMatrix1D       38.0 s          38.0 s             1
TPQ                   1.47 s          1.47 s             1
HardCore             0.179 s         0.179 s             4
ParaH 108            0.184 s         0.110 s             6
ParaH 256            0.489 s         0.392 s             2
ParaH 500            0.404 s         0.327 s             2
ParaH 864            0.324 s         0.267 s             3
Q_TIP4P 1            0.644 s         0.644 s             1
Q_TIP4P 2             1.27 s          1.27 s             1
Q_TIP4P 4             1.28 s          1.28 s             1
HardCore_cuda         36.8 s          36.8 s             1
ParaH cuda 108       0.378 s         0.378 s             2
ParaH cuda 256       0.608 s         0.608 s             1
ParaH cuda 500       0.535 s         0.535 s             1
ParaH cuda 864       0.436 s         0.436 s             2
ParaH auto 108       0.202 s         0.109 s             6
ParaH auto 256       0.295 s         0.238 s             3
ParaH auto 500       0.281 s         0.246 s             3
ParaH auto 864       0.169 s         0.157 s             4
Q_TIP4P cuda         0.991 s         0.991 s             1
Mnist                 14.0 s          14.0 s             1
