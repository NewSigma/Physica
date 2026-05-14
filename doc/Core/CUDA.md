<!--
Copyright 2025-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# CUDA

## Overview

Following C++ design principles, GPU resource management uses RAII. For any host-side object T that requires device parallelism, there exists a template specialization `device_obj<T>` providing:

``` C++
device_obj<T> T::toDevice();
device_obj<T> T::toDeviceAsync();
void T::toDevice(device_obj<T>&);
void T::toDeviceAsync(device_obj<T>&);
```

to copy host data to device;

``` C++
T device_obj<T>::toHost();
T device_obj<T>::toHostAsync();
void device_obj<T>::toHost(T&);
void device_obj<T>::toHostAsync(T&);
```

to copy device data to host.

Only the core logic in the asynchronous non-construction interface needs to be implemented; the other three are wrappers around it:

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
