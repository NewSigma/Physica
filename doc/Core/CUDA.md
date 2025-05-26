# CUDA

## Overview

秉承C++设计原则, 使用RAII进行GPU资源管理。对任意需要device并行的host端对象T, 存在模板特化device_obj<T>使得

``` C++
device_obj<T> T::toDevice();
device_obj<T> T::toDeviceAsync();
void T::toDevice(device_obj<T>&);
void T::toDeviceAsync(device_obj<T>&);
```

将主机数据拷贝到设备;

``` C++
T device_obj<T>::toHost();
T device_obj<T>::toHostAsync();
void device_obj<T>::toHost(T&);
void device_obj<T>::toHostAsync(T&);
```

将设备数据拷贝到主机。

仅需在异步接口非构造中实现核心逻辑, 其余三者在该接口基础上封装:

```C++
auto T::toDevice() const {
    auto result = toDeviceAsync();
    CUDAExecutor::wait();
    return result;
}

auto T::toDeviceAsync() const {
    device_obj<This> result{};
    toDeviceAsync(result);
    return result;
}

void T::toDevice(device_obj<This>& obj) const {
    toDeviceAsync(obj);
    CUDAExecutor::wait();
 }

void T::toDeviceAsync(device_obj<This>& obj) const {
    /* Impl */
}
```
