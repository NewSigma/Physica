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
# TestNVE

This example demonstrates how to use the testNVE function to test the convergence of the integration timestep.

![](./TestNVE.png)

**Figure 1** Using silica as an example: system internal energy as a function of time for different timesteps. Since the NVE ensemble is used, the system energy should be conserved. The 5 fs timestep curve shows energy oscillations due to an excessively large step and is unsuitable for practical simulations. Depending on accuracy requirements, choosing 1 fs or 2 fs is more appropriate.
