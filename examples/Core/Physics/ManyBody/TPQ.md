# TPQ

## Introduction

热力学纯态(Thermo Pure Quantum, TPQ)$^{[1]}$, 用于计算有限温多体系统的各种力学量和热力学量。随着系统格点数增加，计算结果依概率以指数速度收敛到热力学极限。这种收敛是相当快的，实际上4格点的条件下已与精确对角化结果相比无人眼可见差异。TPQ相比精确对角化的改进在于，将计算有限温所需的$O(N)$个本征向量减少到1个随机向量，极大地减少了内存需求。

## 数值结果

![](./TPQ_FD.png)

**图1** 分别使用Physica和HPhi 3.5.2$^{[2]}$计算Hubbard模型倒温度$\beta$与粒子密度间关系，使用全对角化和TPQ两种方法，$U/t = 8$，格点数$L = 4$。由于不能自适应地调整展开阶数，HPhi计算结果不稳定。

    HPhi输入(FullDiag)
    model = "Fermion HubbardGC"
    method = "Full Diag"
    lattice = "chain"
    L = 4
    U = 8
    t = 1

    HPhi输入(TPQ)
    model = "Fermion HubbardGC"
    method = "cTPQ"
    lattice = "chain"
    L = 4
    U = 8
    t = 1
    lanczos_max = 41
    LargeValue = 10
    NumAve = 100
    # ExpandCoef = 10 增加展开阶数可以减小TPQ和FullDiag的偏差

## Reference

[1] Phys. Rev. Lett. 111, 010401 (2013); https://doi.org/10.1103/PhysRevLett.108.240401
[2] HPhi; https://github.com/issp-center-dev/HPhi.git
