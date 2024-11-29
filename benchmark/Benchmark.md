# Benchmark

OS: Ubuntu24.04  
Compiler: clang 17.0.6 + NVCC 12.6.20  
Hardware: Intel(R) Core(TM) i7-14700KF + 32G + NVIDIA GeForce RTX 4060  
Date: 2024-11-29  
Maintainer: Weibo He (NewSigma@163.com)  

----------------------------------------------------------
Benchmark                Time             CPU   Iterations
----------------------------------------------------------
HubbardMatrix1D       5.71 ms         5.71 ms          120
TPQ                   1.11 s          1.11 s             1
lncosh float32         936 ns          936 ns       743537
lncosh float64        2125 ns         2125 ns       327960
lncosh cfloat32       5705 ns         5705 ns       120842
lncosh cfloat64      13446 ns        13446 ns        51884
HardCore              5.19 us         5.19 us       172608
ParaH/108            0.411 ms        0.211 ms         3367
ParaH/256             1.12 ms        0.883 ms          801
ParaH/500             1.82 ms         1.49 ms          474
ParaH/864             2.99 ms         2.50 ms          292
Q_TIP4P               7.01 ms         7.01 ms          100
ParaH cuda/108       0.911 ms        0.911 ms          763
ParaH cuda/256        1.45 ms         1.45 ms          492
ParaH cuda/500        2.53 ms         2.53 ms          276
ParaH cuda/864        4.13 ms         4.13 ms          172
ParaH auto/108       0.422 ms        0.197 ms         3599
ParaH auto/256       0.811 ms        0.631 ms         1126
ParaH auto/500        1.37 ms         1.21 ms          587
ParaH auto/864        2.13 ms         1.97 ms          356
Q_TIP4P cuda          4.09 ms         4.09 ms          171
Mnist                 17.6 s          17.6 s             1
