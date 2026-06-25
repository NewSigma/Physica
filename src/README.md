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
# Project structure

$\quad$The project structure is often the first aspect users encounter when engaging with open-source software. A clean and goal-oriented project organization is as crucial as well-designed and well-documented code. The top-level directory structure of *Physica* follows conventions commonly adopted by open-source C++ projects, consisting of the following directories: `3rdparty`, `benchmark`, `doc`, `examples`, `include`, `src`, and `test`. Within *Physica*, the secondary structure is organized modularly. The modules included in *Physica* are listed in the following table:

| Module | Description |
| -------- | -- |
| Core | Implementation of Physica’s core functionality |
| Gui | Includes 2D and 3D drawing support, using Qt as the drawing backend |
| Logger | A high-performance logging library based on [NanoLog](https://github.com/PlatformLab/NanoLog) |
| Transforms | LLVM Pass Plugin |

$\quad$The secondary project structure extends to each directory within the top-level layout. Tertiary and finer-grained structures encompass APIs and implementation details. The organization follows scientific—rather than purely engineering—logic, allowing domain experts to more readily adopt Physica. Engineering complexities are encapsulated within directories suﬀixed with “Impl”, with lower-level logic nested deeper in the directory hierarchy.

$\quad$Overall, *Physica* employs a goal-oriented, layered structure that progressively exposes complexity. Users not concerned with implementation details will generally encounter fewer such details, as most common use cases are addressed at shallower directory levels. Since scientific workflows are diverse and often require flexibility, we respect users' right to access and modify any underlying details. The pervasive use of templates further facilitates non-intrusive customization and extension.
