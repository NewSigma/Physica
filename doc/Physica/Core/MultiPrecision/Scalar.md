# Scalar的代数结构

定义Scalar的二元运算加和乘(+, *)

## 0的定义

byte最高位为0的任何值(不考虑误差)

## 相反数(-a)的定义

a的相反数(-a)仅与a的length相差一个符号

## 加的代数结构

1.(a + b) + c != a + (b + c)
2.a + 0 = a
3.a + (-a) = 0(待证)
4.a + b = b + a

## 倒数的定义

a ^ -1为a的乘法逆元 (a != 0)

## 乘的代数结构

5.(a * b) * c = a * (b * c)(待证)
6.a * 1 = a
7.a * a ^ -1 = 1 (a != 0)
8.a * b = b * a(待证)
9.a * (b + c) = a * b + a * c(待证)
10.a * b = a * c (a != 0) -> b = c

## 减的定义

a - b = a + (-b)

# About formatted output of Scalar  

std::ostream& operator<<(std::ostream& os, const Scalar<MultiPrecision, false>& s):  

A Scalar s can be expressed in $$s = b \times 2 ^ p \qquad 1 \leq b \le 2 ^ {64} \quad
 -INT\_MAX + 1 \leq p \leq INT\_MAX$$  
 
Express $2^p$ in scientific counting:  

$$2^p = m * 10^n$$  

$$p \ln 2 = \ln m + n \ln 10$$  

let $$n = int(\frac{m \ln 2}{\ln 10})$$  

$$m = e ^ {p \ln 2 - n \ln 10}$$  

Then we have $$s = (b \times m) \times 10^{n}$$
