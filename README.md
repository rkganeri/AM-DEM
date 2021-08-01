# AM-DEM

A discrete element method (DEM) code for simulating powder deposition during additive manufacturing (AM).
A full description of the governing equations and basic methodology is provided in the following paper:

Ganeriwala, R., & Zohdi, T. I. (2016). A coupled discrete element-finite difference model of selective 
laser sintering. Granular Matter, 18(2). https://doi.org/10.1007/s10035-016-0626-0

A prior, open-access paper describing the basics (without the modified Hertzian contact model) is provided 
here:
Ganeriwala, R., & Zohdi, T. I. (2014). Multiphysics modeling and simulation of selective laser sintering 
manufacturing processes. In Procedia CIRP (Vol. 14). https://doi.org/10.1016/j.procir.2014.03.015

This code is based off some sloppy Fortran code I originally wrote in grad school to simulate the process.
The old code is contained in the original\_code folder, which was used to produce the results of the
2016 paper.

The new code within src was re-written in C++ and is capable of running on GPUs via the Kokkos offloading
library.  If Kokkos and this code are compiled with either CUDA or HIP the calculations can be performed
on either Nvidia or AMD GPUs, respectively.  However, if only CPUs are available when compiling, 
Kokkos still provides parallelism via OpenMP threading.  Note that my 2011 MB Pro laptop I've been running/
testing these simulations on does not have any compute capable GPUs so all testing has been CPU-only. 
While I have coded this in a manner such that it should run correctly on GPUs (without unified memory)
as well, I have not tested that so there's a chance some small bugs are present.

Note that Kokkos automatically optimizes memory usage for multi-dimensional arrays based upon the compilation
configuration (e.g. LayoutRight memory for CPU-only or LayoutLeft for GPUs).  Additionally, using Kokkos
we are able to use either OpenMP or CUDA/HIP (depending on compilation setting) without the need for
any #ifdef's within the code!

Feel free to shoot me a message with any bug reports or suggestions for improvement.

