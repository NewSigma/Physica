<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2026 Weibo He
SPDX-License-Identifier: GPL-3.0-or-later
-->

# Physica Core Agent Guide

## Conventions

1. Keep host and device public interfaces in step:
    - For every host type `T` with a `device_obj<T>` specialization, each public member must have a counterpart on the other side: same name, same parameter order and meaning. This includes constructors, operators, getters, static factories, and type aliases.
    - Mirror private helpers, data members, and declaration order where practical, so both implementations can be reviewed side by side.
    - A member may exist on only one side when it is inherently host-only or device-only (for example `__device__`-only entry points or host-side I/O).
