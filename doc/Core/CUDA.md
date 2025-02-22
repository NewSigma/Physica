# CUDA

## Overview

秉承C++设计原则, 使用RAII进行GPU资源管理。对任意需要device并行的host端对象(主机对象)T, 存在模板特化device_obj<T>使得

```
device_obj<T> T::toDevice();
device_obj<T> T::toDeviceAsync();
void T::toDevice(device_obj<T>&);
void T::toDeviceAsync(device_obj<T>&);
```

将主机数据拷贝到设备;

```
T device_obj<T>::toHost();
T device_obj<T>::toHostAsync();
void device_obj<T>::toHost(T&);
void device_obj<T>::toHostAsync(T&);
```

将设备数据拷贝到主机。

为充分利用并行能力, 每个线程至少拥有两个CUDA流

1. 数据流 - Default stream
2. 任务流 - CUDAContext::CUDAStream

上述数据拷贝函数由数据流负责, 所有CUDA核函数将被提交到任务流. 数据流与任务流应进行必要的显式同步。
