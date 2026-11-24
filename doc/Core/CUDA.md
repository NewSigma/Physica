<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2025-2026 Weibo He
SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
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
