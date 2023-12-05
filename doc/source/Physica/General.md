# 通用编程规范

参考[1]中的swap函数原型，进一步规定不可与自身发生交换，即任意对象T的swap函数的一般实现为

void T::swap(T& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    /* Details */
}

实践中发现Self swap常导致Bug，通过添加__restrict帮助编译器优化和提供编译期Self swap警告

## Reference

[1] Scott Meyers 著，侯捷 译. Effective C++：改善程序与设计的55个具体做法[M]. 北京: 电子工业出版社, 2011
