# ParaH

测试项目: ParaH, ParaH_cuda, ParaH_auto

## 平台1

测试环境: Ubuntu20.04 + Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz + 16G + NVIDIA GeForce GTX 1660 Ti Mobile
测试时间: 2023-12-21

$$\begin{matrix}
\mathrm{Platform} && \mathrm{CPU/s} && \mathrm{GPU/s} && \mathrm{CPU + GPU/s} \\
108\ \mathrm{Atom} && 0.47(3) && 0.432(4) && 0.26(2) \\
256\ \mathrm{Atom} && 1.36(5) && 0.82(1) && 0.47(1) \\
500\ \mathrm{Atom} && 1.06(3) && 0.742(4) && 0.39(1) \\
864\ \mathrm{Atom} && 0.86(2) && 0.565(9) && 0.396(4) \\
\end{matrix}$$

## 平台2

测试环境: Ubuntu20.04 + Intel(R) Xeon(R) CPU E3-1231 v3 @ 3.40GHz + 16G + GeForce GTX 960
测试时间: 2023-10-27

$$\begin{matrix}
\mathrm{Platform} && \mathrm{CPU/s} && \mathrm{GPU/s} && \mathrm{CPU + GPU/s} \\
108\ \mathrm{Atom} && 0.41(1) && 0.767(1) && 0.2783(2) \\
256\ \mathrm{Atom} && 1.252(4) && 1.473(7) && 0.631(2) \\
500\ \mathrm{Atom} && 1.04(3) && 1.408(4) && 0.608(1) \\
864\ \mathrm{Atom} && 0.856(4) && 0.967(4) && 0.412(2) \\
\end{matrix}$$
