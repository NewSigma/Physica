<!--
Copyright 2024 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# DenseVector

## 模板参数

**Length**

Length > 0在栈上构造, Length = 0在堆上构造

**Owner**

仅在设备右值向量中出现，用于标识模板表达式在主机端或设备端构造。设备右值向量使用union实现，使用Owner区分使用的成员
