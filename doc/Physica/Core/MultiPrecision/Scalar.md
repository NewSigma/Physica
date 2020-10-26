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
